function [gsymb,sid,tcnt,ccnt,uniq_locs]=load_read_counts(fname,nt,nc)
%function [gsymb,sid,tcnt,ccnt,uniq_locs]=load_read_counts(fname,nt,nc)
%
%RATIONALE: this script imports screen readout in the form of tab
%separated value (TSV) files, with fields: gene_name  shrna_id  treatmentrc controlrc
%Here, gene_name and shrna_id are strings, and *rc are integers.
%Preface a line by # to comment it out.
%
%IN:    fname is the name of the TSV file containing screen readout
%       nt is the number of replicates in treatment
%       nc is the number of replicates in control
%       
%OUT:   
%       gsymb is a cell array of strings containing the gene_name field
%       sid is a cell string array of shrna ids
%       tcnt and ccnt are the counts from the treatment and control files resp. 
%         (uniq_locs==i) is true for the indices of tcnt and ccnt
%         corresponding to gene i (e.g. tcnt(uniq_locs==i) are the shRNA read
%         counts in treatment and sid(uniq_locs==i) are their id strings).

%note: at the moment I am not using inter-replicate variance estimates for
%anything, but this will change in future releases
h=waitbar(0.25,['reading ' regexprep(fname,'_','\\_') '...']);
f=fopen(fname);
s='%s%s';for i=1:(nt+nc), s=[s,'%n']; end
E=textscan(f,s,'Delimiter','\t','CommentStyle','#','HeaderLines',1);
flg=1;for i=1:length(E)-1, flg=flg&&length(E{i})==length(E{i+1});end
if ~flg
    delete(h)
    alert('title','File parse error!','String',['Error parsing line ' num2str(length(E{1}))]);
    return;
end
[gsymb,~,uniq_locs]=unique(E{1});
sid=E{2};
tmp=[];
for i=1:nc
    tmp=[tmp,E{2+i}];
end
if nc>1, ccnt=median(tmp')';
else, ccnt=tmp; end
tmp=[];
for i=1:nt
    tmp=[tmp,E{2+nc+i}];
end
if nt>1, tcnt=median(tmp')';
else, tcnt=tmp; end
waitbar(1,h,'Finishing!')
delete(h)