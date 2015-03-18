function varargout = pareto_rank_tool(varargin)
% PARETO_RANK_TOOL MATLAB code for pareto_rank_tool.fig
%      PARETO_RANK_TOOL, by itself, creates a new PARETO_RANK_TOOL or raises the existing
%      singleton*.
%
%      H = PARETO_RANK_TOOL returns the handle to a new PARETO_RANK_TOOL or the handle to
%      the existing singleton*.
%
%      PARETO_RANK_TOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARETO_RANK_TOOL.M with the given input arguments.
%
%      PARETO_RANK_TOOL('Property','Value',...) creates a new PARETO_RANK_TOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pareto_rank_tool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pareto_rank_tool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pareto_rank_tool

% Last Modified by GUIDE v2.5 18-Mar-2015 11:58:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pareto_rank_tool_OpeningFcn, ...
                   'gui_OutputFcn',  @pareto_rank_tool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pareto_rank_tool is made visible.
function pareto_rank_tool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pareto_rank_tool (see VARARGIN)

% Choose default command line output for pareto_rank_tool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if length(varargin)>0,smp=varargin{1};else,return;end
main_data=get(handles.pareto_rank_root,'UserData');
main_data.algo='pareto';%ranking algorithm
main_data.norm='sd';%normalization scheme
%main_data.pcnt=smp.pcnt;
%main_data.ncnt=smp.ncnt;
%main_data.gsymb=smp.gsymb;
%main_data.uniq_locs=smp.uniq_locs;
main_data.smp=varargin{1};
main_data.sample_list=varargin{2};
set(handles.pareto_rank_root,'UserData',main_data);

% --- Outputs from this function are returned to the command line.
function varargout = pareto_rank_tool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.pareto_rank_root,'UserData');
sample_data=get(main_data.sample_list,'UserData');
n=length(sample_data.gsymb);
%by default normalize by seq depth
pc=sample_data.pcnt;nc=sample_data.ncnt;
sdsf=sum(pc)/sum(nc);
nc=nc*sdsf;
if strcmp(main_data.norm,'uc')
    [fname,pname]=uigetfile('*.*','Select a normalization file...');
    if ~isstr(fname),return; end
    h=waitbar(0,['loading ' regexprep(fname,'_','\\_') ' ...']);
    f=fopen(fullfile(pname,fname));
    %the assumed format is TSV with fields: norm_seq_grp norm_seq_id
    %control_read_count treatment_read_count
    %norm_seq_grp and norm_seq_id are strings and read_count is numeric
    D=textscan(f,'%s%s%n%n','Delimiter','\t','CommentStyle','#','HeaderLines',1);
    fclose(f);
    tmed=median(D{4});
    cmed=median(D{3}*sdsf);
    nc=nc*tmed/cmed;
    waitbar(1,h,'Done!');
    delete(h)
end
h=waitbar(0.25,'Computing ranking...');
slt=get(handles.memory_menu,'Value');
if slt~=3
    clus=parcluster('local');
    nw=clus.NumWorkers;
    if slt==2, parpool(max(1,floor(nw/2)));
    else, parpool; end
    [~,~,rnk,mlodz,nsh,ntar,ndes] = pareto_rank_par(pc,nc,n,sample_data.uniq_locs);
else
    [~,~,rnk,mlodz,nsh,ntar,ndes] = pareto_rank(pc,nc,n,sample_data.uniq_locs);
end
sample_data.prank=rnk;sample_data.mlodz=mlodz;
sample_data.nsh=nsh;sample_data.ntar=ntar;
sample_data.ndes=ndes;
%library swap for FDR estimate
waitbar(0.5,h,'Performing library swap...');
[~,~,nrnk,~,~,~,~] = pareto_rank_par(nc,pc,n,sample_data.uniq_locs);
delete(gcp('nocreate'));
waitbar(0.75,h,'Computing FDR...');
tfdr=zeros(max(rnk),1);
for i=1:length(tfdr)
    tfdr(i)=1/n+length(find(rnk<=rnk(i)&nrnk<=rnk(i)))/length(find(rnk<=rnk(i)));
end
for i=1:length(tfdr),tfdr(i)=min(tfdr(rnk>=rnk(i)));end %compute q-value
sample_data.fdr=ones(size(sample_data.mlodz));
for i=1:length(sample_data.gsymb)
    if ~(sample_data.mlodz(i)<=0||sample_data.nsh(i)<1)
        sample_data.fdr(i)=tfdr(i);
    end
end%annotate genes
waitbar(1,h,'Done!');
%sample_data.sample_name=main_data.sample_name;
set(main_data.sample_list,'UserData',sample_data);
delete(h);
[fname,pname]=uiputfile('screen_results.tsv','Select a file to write to...');
if ~isstr(fname)|isempty(fname),close(handles.pareto_rank_root);end
f=fopen(fullfile(pname,fname),'w');
fprintf(f,'gene\trank\tfdr\teffect_size\t#_active_guide-RNA\n');
[~,sidx]=sort(sample_data.prank);
h=waitbar(0,'Writing results to file...');
for i=1:length(sample_data.gsymb)
    waitbar(i/length(sample_data.gsymb),h,'Writing results to file...');
    fprintf(f,'%s\t',sample_data.gsymb{sidx(i)});
    fprintf(f,'%i\t',sample_data.prank(sidx(i)));
    fprintf(f,'%i\t',sample_data.fdr(sidx(i)));
    fprintf(f,'%g\t',sample_data.mlodz(sidx(i)));
    fprintf(f,'%i\n',sample_data.nsh(sidx(i)));
end
delete(h);
fclose(f);
close(handles.pareto_rank_root);

% --- Executes on selection change in norm_menu.
function norm_menu_Callback(hObject, eventdata, handles)
% hObject    handle to norm_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns norm_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from norm_menu
main_data=get(handles.pareto_rank_root,'UserData');
contents = cellstr(get(hObject,'String'));
s=contents{get(hObject,'Value')};
if strcmp(s,'Sequencing depth (default)'), main_data.norm='sd';
elseif strcmp(s,'User control sequences'), main_data.norm='uc';
end
set(handles.pareto_rank_root,'UserData',main_data);

    
% --- Executes during object creation, after setting all properties.
function norm_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to norm_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'Sequencing depth (default)','User control sequences'});


% --- Executes on selection change in memory_menu.
function memory_menu_Callback(hObject, eventdata, handles)
% hObject    handle to memory_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns memory_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from memory_menu


% --- Executes during object creation, after setting all properties.
function memory_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to memory_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
