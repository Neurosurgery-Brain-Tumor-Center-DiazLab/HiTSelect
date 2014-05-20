function D=report_enrichment(genes,id_type,gidx,hnames,hvals,varargin)
%function D=report_enrichment(genes,id_type,gidx,hnames,hvals,varargin)
%
%IN: genes is a vector of gene entrez ids or a cell array of
%        strings containing mRNA refseq ids or gene symbols to print out
%    id_type is a string, either 'entrez', 'refseq', or 'gsymb'
%    outfile is the file name of the file to write to
%    gidx is a vector of entrez ids if id_type=='entrez' or a cell array of
%       strings if id_type=='refseq' or 'symb'
%    hnames is a cell array of strings containing field header names for
%       fields to print out
%    hvals is a length(gidx)-by-length(hnames) numeric matrix containing
%       field values to print out
%    The following options can be added as additional arguement pairs
%    ('Option',value):
%      'OutFile' : string containing a file name to write the results to
%      'Gsymb'   : a cell array of length(gidx) official gene symbols to display
%      'Entrez'  : a numeric vector of length(gidx) of Entrez Ids to display
%      'RefSeq'  : a cell array of length(gidx) of RefSeq accession Ids to display
%
%OUT: D is a matrix of the vaules which have been written to file

options = containers.Map({'outfile','gsymb','entrez','refseq'},{[],[],[],[]});
%# read the acceptable names
optionNames = options.keys;
%# count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('report_enrichment needs propertyName/propertyValue pairs')
end
for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   inpName = lower(pair{1}); %# make case insensitive
   if any(strmatch(inpName,optionNames))
      options(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end
n=length(genes);t=zeros(n,1);fnd=ones(n,1);
if strmatch(lower(id_type),'entrez') 
    for i=1:n
        tmp=find(gidx==genes(i));
        if ~isempty(tmp)
            genes2{i}=num2str(genes(i));
            t(i)=tmp;
        else
            fnd(i)=0;
        end
    end
    genes=genes2;
else
    for i=1:n
        tmp=find(strcmpi(gidx,genes(i)));
        if ~isempty(tmp), t(i)=tmp; 
        else, fnd(i)=0;end
    end
end
%requested genes
frmt='%s';
D{1,1}='Gene';
n=length(t);fnd=find(fnd);
for i=1:n,D{1+i,1}=gene{fnd(i)};end
col=2;
%gene symbols
if ~isempty(options('gsymb'))
    a=options('gsymb');
    frmt=[frmt '%s'];
    D{1,col}='Gene_Symbol';
    for i=1:n,D{1+i,col}=a{t(i)};end
    col=col+1;
end
%entrez symbols
if ~isempty(options('entrez'))
    a=options('entrez');
    frmt=[frmt '%s'];
    D{1,col}='Entrez_ID';
    for i=1:n,D{1+i,col}=num2str(a(t(i)));end
    col=col+1;
end
%RefSeq Ids
if ~isempty(options('refseq'))
    a=options('refseq');
    frmt=[frmt '%s'];
    D{1,col}='RefSeq_ID';
    for i=1:n,D{1+i,col}=a{i};end
    col=col+1;
end
for i=1:length(hnames)
    frmt=[frmt '%s'];
    D{1,col}=hnames{i};
    for j=1:n,D{1+j,col}=num2str(hvals(j,i));end
    col=col+1;
end
frmt=[frmt '\n'];
if ~isempty(option('outfile'))
    f=fopen(option('outfile'),'w');
    for i=1:size(D,1)
        fprintf(f,frmt,D{i,:});
    end
end






