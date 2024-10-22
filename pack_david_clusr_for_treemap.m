function x=pack_david_clusr_for_treemap(c,smp,GO)
%function x=pack_david_clusr_for_treemap(c,smp,GO)
%
%IN:c is an object containing data from a DAVID gene ontology term cluster
%report, generated by query_david.m
%smp is an object with arrays of data concerning the genes
%    see pareto_gui.m
%GO is a gene ontology mapping generated by geneont from a gene ontology
%obo file
%
%OUT:
%   x is a structure encoding the cluster report, packed into a format 
%readable by the InvoVis toolbox
if isdeployed
    f=fopen(fullfile(ctfroot,'david_cluster_report.txt'),'w');
else
    f=fopen(fullfile(pwd,'david_cluster_report.txt'),'w');
end
x=struct('id','root','name','root_node','data',[]); %root node
%collect all the scores and cluster sizes
for i=1:length(c)
    scz(i)=c(i).getScore;
    rec=c(i).getSimpleChartRecords;
    ntm(i)=length(rec);
end
mscz=max(scz);mntm=max(ntm);
%color the GO term clusters by score, size them by cluster size
%a range of colors from blue to red, one for each decile of score/p-value
cmp={'#0000FF','#1900E6','#3300CC','#4D00B2','#800080','#990066','#B2004C','#CC0033','#E6001A','#FF0000'};
%create a node for each cluster (clus) and set them to be the children of x
for i=1:min(length(c),10)
   if exist('d','var'),clear d;end
   rec=c(i).getSimpleChartRecords;%get all the terms in the cluster
   %set the nodes size by the number of GO terms, as a % of the max
   d.area=max(1,round(ntm(i)/mntm*10));%1<=dim,area<=10
   %set the node color as a hex 
   if scz(i) <=.2, pr=1;
   elseif scz(i)>.2&&scz(i)<=.4, pr=2;
   elseif scz(i)>.4&&scz(i)<=.6, pr=3;
   elseif scz(i)>.6&&scz(i)<=.8, pr=4;
   elseif scz(i)>.8&&scz(i)<=1, pr=5;
   elseif scz(i)>1&&scz(i)<=5, pr=6;
   elseif scz(i)>5&&scz(i)<=10, pr=7;
   elseif scz(i)>10&&scz(i)<=15, pr=8;
   elseif scz(i)>15&&scz(i)<=20, pr=9;
   elseif scz(i)>20, pr=10; end
   d.color = cmp{pr};
   d.score = scz(i);
   %set the cluster name to be the most frequently annotated term
   szs=zeros(length(rec),1);ezs=szs;
   for j=1:length(rec),szs(j)=rec(j).getPercent;ezs(j)=rec(j).getEase;end
   [~,tsidx]=max(szs)
   if length(tsidx)==1, tidx=tsidx;
   else
       [~,tidx]=max(ezs(tsidx))
       tidx=min(tsidx(tidx))
   end
   nm=char(rec(tidx).getTermName);
   if nm(1)=='G' %if the cluster name is a GO term, deconstruct
       tstr=textscan(nm,'GO:%s%s','delimiter','~');
       nm=tstr{2}{:};
       d.go=str2num(tstr{1}{:});
       id=['clus_' num2str(i) '_GO:' tstr{1}{:}];
   else
       id=['clus_' num2str(i) '_' nm];
   end
   if ~isempty(GO) %if we have GO information and the cluster name is a
       if isfield(d,'go') %GO term, annotate the node with GO term info
           try
               trm=GO(d.go).terms(1);
               if ~isempty(trm.ontology)
                tmp=strrep(trm.ontology,'"','');
                tmp=strrep(tmp,':','');
                d.ont=tmp;
               end
               if ~isempty(trm.definition)
                tmp=strrep(trm.definition,'"','');
                tmp=strrep(tmp,':','');
                d.def=tmp;
               end
           catch me
           end
       end
   end
   clus=struct('id',id,'name',nm,'data',d);
   clear d;
   %create a node for each GO term in cluster i
   %gtrm are the children of clus(i), one for each GO term
   for j=1:length(rec),pvl(j)=rec(j).getEase;end,mpvl=max(1-pvl);
   for j=1:length(rec)
       if exist('d','var'),clear d;end
       %set the size of the node by the number of genes
       d.area=max(1,round(rec(j).getPercent));
       %set the color by p-value
       if pvl(i) >=.4, pr=1;
       elseif pvl(i)<.4&&pvl(i)>=.3, pr=2;
       elseif pvl(i)<.3&&pvl(i)>=.2, pr=3;
       elseif pvl(i)<.2&&pvl(i)>=0.1, pr=4;
       elseif pvl(i)<0.1&&pvl(i)>=0.05, pr=5;
       elseif pvl(i)<0.05&&pvl(i)>=1e-3, pr=6;
       elseif pvl(i)<1e-3&&pvl(i)>=1e-5, pr=7;
       elseif pvl(i)<1e-5&&pvl(i)>=1e-7, pr=8;
       elseif pvl(i)<1e-7&&pvl(i)>=1e-9, pr=9;
       elseif pvl(i)<1e-9, pr=10; end
       d.color = cmp{pr};
       nm=char(rec(j).getTermName);
       d.gns=char(rec(j).getGeneIds);
       D=textscan(d.gns,'%n','Delimiter',',');
       fprintf(f,'"%s"\t',strrep(nm,',',':'));fprintf(f,'%g\t',pvl(j));fprintf(f,'%g\t',double(rec(j).getPercent));
       for k=1:length(D{1})
           tmpg=find(D{1}(k)==smp.gid);
           if isempty(tmpg), continue, end
           fprintf(f,'%s,',smp.gsymb{tmpg});
       end
       fprintf(f,'\n');
       if strcmp(nm(1:3),'GO:') %if the cluster name is a GO term, deconstruct
           tstr=textscan(nm,'GO:%s%s','delimiter','~');
           nm=tstr{2}{:};
           d.go=str2num(tstr{1}{:});
           id=['trm_' num2str(i) '_' num2str(j) '_GO:' tstr{1}{:}];
       else
           id=['trm_' num2str(i) '_' num2str(j) '_' nm];
       end
       if ~isempty(GO) %if we have GO information and the term name is a
           if isfield(d,'go') %GO term, annotate the node with GO term info
             try
               trm=GO(d.go).terms(1);
               if ~isempty(trm.ontology)
                   tmp=strrep(trm.ontology,'"','');
                   tmp=strrep(tmp,':','');
                   d.ont=tmp;
               end
               if ~isempty(trm.definition)
                   tmp=strrep(trm.definition,'"','');
                   tmp=strrep(tmp,':','');
                   d.def=tmp;
               end
             catch me
             end
           end
       end
       d.gns=char(rec(j).getGeneIds);
       d.pval=pvl(j);
       gtrm=struct('id',id,'name',nm,'data',d);
       %parse the gene ids
       gns=textscan(char(rec(j).getGeneIds),'%s','EndOfLine',',');gns=gns{1};
       d.area=d.area/length(gns);
       for k=1:length(gns),idx(k)=min(find(smp.gid==str2num(gns{k})));end
       if isfield(smp,'prank')&&~isempty(smp.prank),mxr=max(smp.prank);end
       %create a child node of gtrm
       for k=1:length(gns)
           %extract the screen enrichment data and save in the gene nodes
           gstr=smp.gsymb{idx(k)};
           if isfield(smp,'fdr')&&~isempty(smp.fdr)
               d.rank=smp.prank(idx(k));
               d.pvl=smp.fdr(idx(k));
               d.nsh=smp.nsh(idx(k));
               d.mlodz=smp.mlodz(idx(k));
               pct=d.rank/max(smp.prank(idx));
               if pct >=.7, pr=1;
               elseif pct<.7&&pct>=.6, pr=2;
               elseif pct<.6&&pct>=.5, pr=3;
               elseif pct<.5&&pct>=0.4, pr=4;
               elseif pct<0.4&&pct>=0.3, pr=5;
               elseif pct<0.3&&pct>=0.2, pr=6;
               elseif pct<0.2&&pct>=0.1, pr=7;
               elseif pct<0.1&&pct>=0.05, pr=8;
               elseif pct<0.05&&pct>=0.01, pr=9;
               elseif pct<0.01, pr=10; end
               d.color=cmp{pr};
           end
           if isfield(smp,'fc')&&~isempty(smp.fc)
               d.fc=smp.fc(idx(k));
               d.tt=smp.tt(idx(k));
           end
           gtrm.children(k)=struct('id',...
                   ['gn_' num2str(i) '_' num2str(j) '_' num2str(k) '_' gns{k}],...
                   'name',gstr,'data',d);
       end
       clus.children(j)=gtrm;
   end
   x.children(i)=clus;
end
fclose(f);