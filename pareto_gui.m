function varargout = pareto_gui(varargin)
% PARETO_GUI MATLAB code for pareto_gui.fig
%      PARETO_GUI, by itself, creates a new PARETO_GUI or raises the existing
%      singleton*.
%
%      H = PARETO_GUI returns the handle to a new PARETO_GUI or the handle to
%      the existing singleton*.
%
%      PARETO_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARETO_GUI.M with the given input arguments.
%
%      PARETO_GUI('Property','Value',...) creates a new PARETO_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pareto_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pareto_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pareto_gui

% Last Modified by GUIDE v2.5 31-Oct-2013 12:54:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pareto_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @pareto_gui_OutputFcn, ...
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


% --- Executes just before pareto_gui is made visible.
function pareto_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pareto_gui (see VARARGIN)

%splash('*.png')

% Choose default command line output for pareto_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
main_data=get(handles.root_window,'UserData');
main_data.java_loaded=0;
main_data.GO=[];
set(handles.root_window,'UserData',main_data);
%M=imread('pareto.jpg','jpg');
%image(M,'Parent',handles.axes1);
%set(handles.axes1,'ytick',[],'xtick',[]);

% UIWAIT makes pareto_gui wait for user response (see UIRESUME)
% uiwait(handles.root_window);



% --- Outputs from this function are returned to the command line.
function varargout = pareto_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in restore_session_pushbutton.
function restore_session_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to restore_session_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sample_data=get(handles.sample_list,'UserData');
sample_data.main_glst=handles.working_gene_listbox;
set(handles.sample_list,'UserData',sample_data);
main_data=get(handles.root_window,'UserData');
if ~isfield(main_data,'last_dir')
    [fname pname]=uigetfile('*.mat','Restore session from file...');
else
    [fname pname]=uigetfile([main_data.last_dir '*.mat'],'Restore session from file...');
end
if ~isstr(fname),return;end
if ~exist('sample_data','var')
    alert('title','No samples found!','string',['No samples found in ' fname])
else
    load([pname fname],'sample_data');
    set(handles.sample_list,'String',sample_data.smp_lst,'UserData',sample_data);
    alert('title','Sample loaded!','string',['Sample ' fname ' loaded!'])
end
if isfield(sample_data,'glists')
    set(handles.gene_lists,'UserData',sample_data.glists,'String',sample_data.glist_ids);
end
main_data.last_dir=pname;
set(handles.root_window,'UserData',main_data);
set(handles.sample_list,'Value',1);
%set(handles.sample_genome_textbox,'String',sample_data.genome{1});


% --- Executes during object creation, after setting all properties.
function main_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to main_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String',{'Welcome!'});

% --- Executes on button press in export_session_pushbutton.
function export_session_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to export_session_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.root_window,'UserData');
sample_data=get(handles.sample_list,'UserData');
glist_ids=get(handles.gene_lists,'String');
glists=get(handles.gene_lists,'UserData');
if ~isfield(main_data,'last_dir')
    [fname pname]=uiputfile('*.mat','Save data as...','Untitled.mat');
else
    [fname pname]=uiputfile([main_data.last_dir '*.mat'],'Save data as...','Untitled.mat');
end
if ~pname,return;end
main_data.last_dir=pname;
set(handles.root_window,'UserData',main_data);
sample_data.glist_ids=glist_ids;sample_data.glists=glists;
save([pname fname],'sample_data');
alert('title','Working sample saved','string',['Wrote matlab binary file ' fname]);


% --- Executes on button press in ss_sample_pushbutton.
function load_sample_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to load_sample_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.root_window,'UserData');
sample_data=get(handles.sample_list,'UserData');   
sample_data.main_glst=handles.working_gene_listbox;
set(handles.sample_list,'UserData',sample_data);
smp_info=choose_file_type('title','Load data','String','Import data...');
if isempty(smp_info)||isempty(smp_info.sample_type) return; end
fname=[];
if strcmp(smp_info.sample_type,'Screen data')
    if ~isfield(smp_info,'num_treat')||isempty(smp_info.num_treat), return; end
    main_data.num_treat=smp_info.num_treat;
    main_data.num_ctr=smp_info.num_ctr;
    if ~isfield(main_data,'last_dir')
        [fname, pname]=uigetfile('*.*','Select the data file to import...');
    else
        [fname, pname]=uigetfile(fullfile(main_data.last_dir,'*.*'),'Select the data file to import...');
    end
    if ~isstr(fname),return;
    else
        main_data.last_dir=pname;
        set(handles.root_window,'UserData',main_data);
    end
    try
        [gsymb,sid,pcnt,ncnt,uniq_locs]=load_read_counts(fullfile(pname,fname),main_data.num_treat,main_data.num_ctr);
    catch me
        alert('title','Error parsing file...','String','Make sure your file matches specified format.');
        keyboard()
    end
    if strcmp(smp_info.build,'hg19'),load('symb2entrez_hsa.mat');symb2ent=symb2ent_hsa;
    else,load('symb2entrez_mmu.mat');symb2ent=symb2ent_mmu;end
    sample_data.gsymb=gsymb;sample_data.sid=sid;sample_data.pcnt=pcnt;
    sample_data.ncnt=ncnt;sample_data.uniq_locs=uniq_locs;
    sample_data.build=smp_info.build;
    kz=symb2ent.keys;
    sample_data.gid=zeros(size(gsymb));
    r=1;not_found={};
    for i=1:length(gsymb) %look for each gene symbol's entrez id
        idx=find(strcmpi(gsymb{i},kz));
        if ~isempty(idx)
            sample_data.gid(i)=symb2ent(kz{min(idx)});
        else
            not_found{r}=gsymb{i};r=r+1;
        end                         
    end
    sample_data.gexp_pcut=0.05;
    try
        smp_id=set_sample_id('title','Enter sample ID:','string','Enter a name for the sample');
        sample_data.sample_name=smp_id;
        smp_id=['screen data: ',smp_id];
    catch me
        smp_id='screen data';
        sample_data.sample_name=smp_id;
    end
    sample_data.prank=[];sample_data.fdr=[];
    sample_data.smp_lst=smp_id;
elseif strcmp(smp_info.sample_type,'Gene expression data')
    if isempty(sample_data.gsymb),alert('title','No screen data loaded','string','Please load screen data first');end
    if ~isfield(main_data,'last_dir')
        [fname pname]=uigetfile('*.*','Select gene expresion data file...');
    else
        [fname pname]=uigetfile(fullfile(main_data.last_dir,'*.*'),'Select the gene expression data file...');
    end
    if ~isstr(fname),return;
    else
        main_data.last_dir=pname;
        set(handles.root_window,'UserData',main_data);
    end
    f=fopen(fullfile(pname,fname));
    %the assumed format is TSV with fields: gene_name fold_change p-value
    %gene_name is a string, fold_change and p-value are doubles
    D=textscan(f,'%s%n%n');
    h=waitbar(0,'annotating your screen data...');
    for i=1:length(D{1})
        waitbar(i/length(D{1}),h,'annotating your screen data...');
        idx=find(strcmpi(sample_data.gsymb,D{1}{i}));
        if ~isempty(idx)
            sample_data.fc(idx)=D{2}(i);
            sample_data.tt(idx)=D{3}(i);
        end
    end
    delete(h)
    try
        smp_id=set_sample_id('title','Enter sample ID:','string','Enter a name for the sample');
        smp_id=['gene expression data: ',smp_id];
    catch me
        smp_id='gene expression data';
    end
    sample_data.smp_lst={sample_data.smp_lst,smp_id};
end
set(handles.sample_list,'String',sample_data.smp_lst,'UserData',sample_data);
set(handles.root_window,'UserData',main_data);

% --- Executes during object creation, after setting all properties.
function sample_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
sample_data.gsymb=[]; sample_data.sid=[];
sample_data.seq=[]; sample_data.pcnt=[];
sample_data.ncnt=[]; sample_data.uniq=[];
sample_data.uniq_locs=[];sample_data.n=[];
sample_data.prank=[];sample_data.fname=[];
sample_data.nsh=[];sample_data.mlodz=[];
sample_data.fc=[];sample_data.tt=[];
sample_data.full_path={};sample_data.smp_lst={};
sample_data.type={};sample_data.genome={};
sample_data.net_cent_comb=[];
set(hObject,'UserData',sample_data);

% --- Executes on button press in delete_sample_pushbutton.
function delete_sample_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to delete_sample_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sample_data=get(handles.sample_list,'UserData');
rmidx=get(handles.sample_list,'Value');
%contents = cellstr(get(handles.sample_list,'String'));
%to_rm=contents{rmidx};
if rmidx==1
    set(handles.sample_list,'String',{});
    sample_data.gsymb=[]; sample_data.sid=[];
    sample_data.seq=[]; sample_data.pcnt=[];
    sample_data.ncnt=[]; sample_data.uniq=[];
    sample_data.uniq_locs=[];sample_data.n=[];
    sample_data.prank=[];sample_data.fname=[];
    sample_data.nsh=[];sample_data.mlodz=[];
    sample_data.fc=[];sample_data.tt=[];
    sample_data.type={};sample_data.genome={};
    sample_data.full_path={};sample_data.smp_lst='';
    set(hObject,'UserData',sample_data);
elseif rmidx==2
    sample_data.fc=[];sample_data.tt=[];
    sample_data.full_path=sample_data.full_path(1);
    sample_data.smp_lst=sample_data.smp_lst{1};
end
set(handles.sample_list,'String',sample_data.smp_lst);
set(hObject,'UserData',sample_data);


% --- Executes on button press in pareto_rank.
function pareto_rank_Callback(hObject, eventdata, handles)
% hObject    handle to pareto_rank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.root_window,'UserData');
sample_data=get(handles.sample_list,'UserData');
if isempty(sample_data.pcnt)
    alert('title','Screen data not loaded','string','Load screen data first...');
    return;
end
pareto_rank_tool(sample_data,handles.sample_list);


% --- Executes on button press in pareto_plot.
function pareto_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pareto_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sample_data=get(handles.sample_list,'UserData');
sample_data.main_glst=handles.working_gene_listbox;
sample_data.glists_textbox=handles.gene_lists;
if ~isfield(sample_data,'prank')||isempty(sample_data.prank)
    alert('String','Perform a gene ranking or restore one from a saved session first')
    return;
end
pareto_plt_gui(sample_data);

% --- Executes on selection change in working_gene_listbox.
function working_gene_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to working_gene_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns working_gene_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from working_gene_listbox
gnz=get(hObject,'UserData');
if isempty(gnz), return; end
s=gnz{get(hObject,'Value')};
sample_data=get(handles.sample_list,'UserData');
glst_ids=get(handles.working_gene_listbox,'UserData');
idx=min(find(strcmpi(s,sample_data.gsymb)));
set(handles.screen_rank_textbox,'String',num2str(sample_data.prank(idx)));
if sample_data.fdr(idx)<=0.05
    set(handles.screen_fdr_textbox,'String',num2str(sample_data.fdr(idx)),'ForegroundColor','r');
else
    set(handles.screen_fdr_textbox,'String',num2str(sample_data.fdr(idx)),'ForegroundColor','k');
end
set(handles.nsh_textbox,'String',num2str(sample_data.nsh(idx)));
if isempty(sample_data.tt)|sample_data.tt(idx)==-1,
    set(handles.ttest_textbox,'String','NA');
    set(handles.exp_fc_textbox,'String','NA');    
else
    if sample_data.tt(idx)<=0.05
        set(handles.ttest_textbox,'String',num2str(sample_data.tt(idx)),'ForegroundColor','r');
    else
        set(handles.ttest_textbox,'String',num2str(sample_data.tt(idx)),'ForegroundColor','k');
    end
    set(handles.exp_fc_textbox,'String',num2str(sample_data.fc(idx)));
end
set(handles.mlodz_textbox,'String',num2str(sample_data.mlodz(idx)));
if isfield(sample_data,'net_cent_comb')&~isempty(sample_data.net_cent_comb)
    set(handles.net_cent_textbox,'String',num2str(sample_data.net_cent_comb(idx)));
else, set(handles.net_cent_textbox,'String','NA');
end

% --- Executes during object creation, after setting all properties.
function working_gene_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to working_gene_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{});
set(hObject,'UserData',{});


% --- Executes on button press in delete_gene_pushbutton.
function delete_gene_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to delete_gene_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx=get(handles.working_gene_listbox,'Value');
glist_old=get(handles.working_gene_listbox,'String');
glist_ids_old=get(handles.working_gene_listbox,'UserData');
if isempty(idx)
    set(handles.working_gene_listbox,'String',{},'UserData',{},'Value',1);
else
    j=1;glist={};glist_ids={};
    for i=1:idx-1,glist{j}=glist_old{i};glist_ids{j}=glist_ids_old{i};j=j+1;end
    for i=idx+1:length(glist_old),glist{j}=glist_old{i};glist_ids{j}=glist_ids_old{i};j=j+1;end
    set(handles.working_gene_listbox,'String',glist,'UserData',glist_ids,'Value',max(1,idx-1));
end
    
function find_gene_editbox_Callback(hObject, eventdata, handles)
% hObject    handle to find_gene_editbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of find_gene_editbox as text
%        str2double(get(hObject,'String')) returns contents of find_gene_editbox as a double

% --- Executes during object creation, after setting all properties.
function find_gene_editbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to find_gene_editbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in find_gene_button.
function find_gene_button_Callback(hObject, eventdata, handles)
% hObject    handle to find_gene_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s=get(handles.find_gene_editbox,'String');
sample_data=get(handles.sample_list,'UserData');
glst_ids=get(handles.working_gene_listbox,'UserData');
idx=find(strcmpi(s,glst_ids));
if isempty(idx),set(handles.find_gene_editbox,'String','Gene not found!');return;end
pidx=min(idx);
idx=min(find(strcmpi(s,sample_data.gsymb)));
if ~isempty(pidx)
    set(handles.working_gene_listbox,'Value',pidx);
    set(handles.screen_rank_textbox,'String',num2str(sample_data.prank(idx)));
    if sample_data.fdr(idx)<=0.05
        set(handles.screen_fdr_textbox,'String',num2str(sample_data.fdr(idx)),'ForegroundColor','r');
    else
        set(handles.screen_fdr_textbox,'String',num2str(sample_data.fdr(idx)),'ForegroundColor','k');
    end
    set(handles.nsh_textbox,'String',num2str(sample_data.nsh(idx)));
    if isempty(sample_data.tt)|sample_data.tt(idx)==-1,
        set(handles.ttest_textbox,'String','NA');
        set(handles.exp_fc_textbox,'String','NA');    
    else
        if sample_data.tt(idx)<=0.05
            set(handles.ttest_textbox,'String',num2str(sample_data.tt(idx)),'ForegroundColor','r');
        else
            set(handles.ttest_textbox,'String',num2str(sample_data.tt(idx)),'ForegroundColor','k');
        end
        set(handles.exp_fc_textbox,'String',num2str(sample_data.fc(idx)));
    end
    set(handles.mlodz_textbox,'String',num2str(sample_data.mlodz(idx)));
    if isfield(sample_data,'net_cent_comb')
        set(handles.net_cent_textbox,'String',num2str(sample_data.net_cent_comb(idx)));
    else
        set(handles.net_cent_textbox,'String','NA');
    end
end



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over working_gene_listbox.
function working_gene_listbox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to working_gene_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in david_analysis_pushbutton.
function david_analysis_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to david_analysis_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sample_data=get(handles.sample_list,'UserData');
main_data=get(handles.root_window,'UserData');
gene_lst=get(handles.working_gene_listbox,'String');
gene_ids=get(handles.working_gene_listbox,'UserData');
if isempty(gene_ids)
    alert('title','Empty working gene list...','String','Please populate the working gene list first (e.g. via Plot screen readout)...');
    return;
end
for i=1:length(gene_ids)
    idx(i)=min(find(strcmpi(gene_ids{i},sample_data.gsymb)));
end
idx=intersect(idx,find(sample_data.gid~=0));
if isempty(sample_data.fc)
    smp.fc=[];
    smp.tt=[];
else
    smp.fc=sample_data.fc(idx);
    smp.tt=sample_data.tt(idx);
end
smp.gid=sample_data.gid(idx);
smp.gexp_pcut=sample_data.gexp_pcut;
smp.fdr=sample_data.fdr(idx);
smp.nsh=sample_data.nsh(idx);
smp.mlodz=sample_data.mlodz(idx);
smp.prank=sample_data.prank(idx);
smp.gsymb=sample_data.gsymb(idx);
smp.pname=main_data.last_dir;
smp.build=sample_data.build;
smp.pareto_gui_root_handle=handles.root_window;
smp.glists_textbox=handles.gene_lists;
smp.main_glst=handles.working_gene_listbox;
david_tool(smp)


% --- Executes on button press in genemania_pushbutton.
function genemania_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to genemania_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sample_data=get(handles.sample_list,'UserData');
main_data=get(handles.root_window,'UserData');
if isempty(sample_data.gsymb),alert('Title','No data loaded','String','Load some data first :)');end
gene_lst=get(handles.working_gene_listbox,'String');
gene_ids=get(handles.working_gene_listbox,'UserData');
if isempty(gene_ids)
    alert('Title','No genes selected','String','Add genes to the Working gene list first');
    return;
end
for i=1:length(gene_ids)
    idx(i)=min(find(strcmpi(gene_ids{i},sample_data.gsymb)));
end
if isempty(sample_data.fc)
    smp.fc=ones(size(idx'));
    smp.tt=ones(size(idx'));
else
    smp.fc=sample_data.fc(idx);
    smp.tt=sample_data.tt(idx);
end
smp.fdr=sample_data.fdr(idx);
smp.nsh=sample_data.nsh(idx);
smp.mlodz=sample_data.mlodz(idx);
smp.prank=sample_data.prank(idx);
smp.gsymb=sample_data.gsymb(idx);
smp.pname=main_data.last_dir;
smp.build=sample_data.build;
smp.gid=sample_data.gid(idx);
smp.pareto_gui_root=handles.root_window;
smp.main_sample_list=handles.sample_list;
genemania_tool(smp)

% --- Executes on selection change in gene_lists.
function gene_lists_Callback(hObject, eventdata, handles)
% hObject    handle to gene_lists (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns gene_lists contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gene_lists
glists_data=get(hObject,'UserData');%get gene lists data
glists_names=get(hObject,'String');%get the name of all the gene lists
idx=get(hObject,'Value');
if ~isempty(idx)
    set(handles.working_gene_listbox,'String',glists_data{idx},'UserData',glists_data{idx},'Value',1);
else
    set(handles.working_gene_listbox,'String',{},'UserData',{},'Value',1);
end



% --- Executes during object creation, after setting all properties.
function gene_lists_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gene_lists (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in new_list_pushbutton.
function new_list_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to new_list_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

glists_data=get(handles.gene_lists,'UserData');%get gene lists data
glists_names=get(handles.gene_lists,'String');%get the name of all the gene lists
glists_data{end+1}={};
lst_name=set_sample_id('title','Enter gene list ID:','string',sprintf(['Enter a name for the gene list.']));
if isempty(lst_name),glists_names{end+1}=['new_list' num2str(length(glists_names))];
else,glists_names{end+1}=lst_name;end
set(handles.gene_lists,'UserData',glists_data,'String',glists_names,'Value',length(glists_names));
set(handles.working_gene_listbox,'UserData',{},'String',{});

% --- Executes on button press in delete_list_pushbutton.
function delete_list_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to delete_list_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
glists_data=get(handles.gene_lists,'UserData');%get gene lists data
glists_names=get(handles.gene_lists,'String');%get the name of all the gene lists
idx=get(handles.gene_lists,'Value')
if length(glists_names)==1,gl_data={};gl_names={};
else
    j=1;
    for i=1:idx-1
        gl_data{j}=glists_data{i};
        gl_names{j}=glists_names{i};
        j=j+1;
    end
    for i=idx+1:length(glists_data)
        gl_names{j}=glists_names{i};
        gl_data{j}=glists_data{i}
        j=j+1;
    end
end
set(handles.gene_lists,'UserData',gl_data,'String',gl_names,'Value',max(1,idx-1));

% --- Executes on button press in export_gene_list_pushbutton.
function export_gene_list_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to export_gene_list_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sample_data=get(handles.sample_list,'UserData');%get screen data
main_data=get(handles.root_window,'UserData');%get global data
glists_data=get(handles.gene_lists,'UserData');%get gene lists data
glists_names=get(handles.gene_lists,'String');%get the name of all the gene lists
if isempty(glists_data),return;end
idx=get(handles.gene_lists,'Value');
gnz=glists_data{idx};
if ~isfield(main_data,'last_dir')
    [fname pname]=uiputfile(glists_names{idx},'Select a file to write to...');
else
    main_data.last_dir
    [fname pname]=uiputfile('*.*','Select a file to write to...',fullfile(main_data.last_dir,[strrep(strrep(glists_names{idx},' ','_'),',',''),'.txt']));
end
if ~isstr(fname),return;
else
    main_data.last_dir=pname;
    set(handles.root_window,'UserData',main_data);
end
f=fopen(fullfile(pname,fname),'w');
fprintf(f,'gene\tscreen_rank\tscreen_fdr\tcollective_hairpin_activity\tnumber_of_active_hairpins');
if isfield(sample_data,'fc')&&~isempty(sample_data.fc)
    fprintf(f,'\tgexp_fold_change\tttest_pvalue');
end
if isfield(sample_data,'net_cent_comb')&&~isempty(sample_data.net_cent_comb)
    fprintf(f,'\tnetwork_centrality');
end
fprintf(f,'\n');
for i=1:length(gnz)
    pidx=min(find(strcmpi(gnz{i},sample_data.gsymb)));
    fprintf(f,'%s\t',gnz{i});
    fprintf(f,'%u\t',sample_data.prank(pidx));
    fprintf(f,'%g\t',sample_data.fdr(pidx));
    fprintf(f,'%g\t',sample_data.mlodz(pidx));
    fprintf(f,'%u',sample_data.nsh(pidx));
    if isfield(sample_data,'fc')&&~isempty(sample_data.fc)
        fprintf(f,'\t%g\t',sample_data.fc(pidx));
        fprintf(f,'%g',sample_data.tt(pidx));
    end
    if isfield(sample_data,'net_cent_comb')&&~isempty(sample_data.net_cent_comb),fprintf(f,'\t%g',sample_data.net_cent_comb(pidx));end
    fprintf(f,'\n');
end



% --- Executes on button press in import_gene_list_pushbutton.
function import_gene_list_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to import_gene_list_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sample_data=get(handles.sample_list,'UserData');%get screen data
if isempty(sample_data.gsymb),alert('String','Load screen data first...'),return;end
main_data=get(handles.root_window,'UserData');%get global data
glists_data=get(handles.gene_lists,'UserData');%get gene lists data
glists_names=get(handles.gene_lists,'String');%get the name of all the gene lists
idx=get(handles.gene_lists,'Value');
%get the filename
if ~isfield(main_data,'last_dir')
    [fname pname]=uigetfile('*.*','Select the gene list...');
else
    [fname pname]=uigetfile([main_data.last_dir,'*.*'],'Select the gene list...');
end
if ~isstr(fname),return;
else
    main_data.last_dir=pname;
    set(handles.root_window,'UserData',main_data);
end
f=fopen(fullfile(pname,fname));
D=textscan(f,'%s',1);
if strcmpi(D{1},'gene'),hl=1;else,hl=0;end%you can have one headerline starting with "gene"
fseek(f,0,-1);%rewind file
D=textscan(f,'%s%*[^\n]','HeaderLines',hl);
glist={};j=1;
for i=1:length(D{1})
    if any(strcmpi(D{1}{i},sample_data.gsymb)), glist{j}=D{1}{i};j=j+1;end
end
set(handles.working_gene_listbox,'String',glist,'UserData',glist,'Value',1);
lst_name=[];
lst_name=set_sample_id('title','Enter gene list ID:','string',sprintf(['Enter a name for the gene list\n(' fname ')']));
if isempty(lst_name),glists_names{end+1}=['new_list' num2str(length(glists_names))];
else,glists_names{end+1}=lst_name;end
glists_data{end+1}=glist;
set(handles.gene_lists,'UserData',glists_data,'String',glists_names,'Value',length(glists_names));




% --- Executes on key press with focus on gene_lists and none of its controls.
function gene_lists_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to gene_lists (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in sample_list.
function sample_list_Callback(hObject, eventdata, handles)
% hObject    handle to sample_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sample_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sample_list
