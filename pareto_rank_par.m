function [M1, M2, rnk, mlodz, nsh, ntar, ndes] = pareto_rank_par(pos_cnt,neg_cnt,n,uniq_locs)
%function [M1 M2 rnk mlodz nsh ntar ndes] =pareto_rank_par(pos_cnt,neg_cnt,n,uniq_locs)
%
%IN: *cnt: shrna read count vectors
%    n is the number of genes in the screen
%    uniq_locs: a vector of indices into unique(*cnt), ie. find(uniq_locs==i) gives the indices of all shrna targeting gene i
%           uniq and uniq_locs are produced by [~,uniq,uniq_locs]=unique(gene_id_list)
%    graph_on is set to 1 if a plot is requested
%
%OUT:
%rnk(i) is the sum of all the probabilities that gene i dominates gene j, summed over j
%plus half the sum of all the probabilities that it is not comparible to
%gene j, summed over j
%mlodz is the genewise weighted mean log odds ratio of ips+ over ips-
%  the wieghts are chosen via a DerSimonian-Laird random effects model
%nsh is the genewise number of active targeting shRNA
%ntar is the genewise number of sequenced targeting shRNA
%ndes is the genewise number of infected (# in original library) targeting shRNA

%add a prior and compute shrna read count odds
prior=1;pc=pos_cnt+prior;nc=neg_cnt+prior;
psum=sum(pc);nsum=sum(nc);
odzp=(pc/psum)./(1-pc/psum);odzn=(nc/nsum)./(1-nc/nsum);
mlodz=zeros(n,1);%genewise mean log odds ratios
ntar=zeros(n,1);%number of sequenced shrna targeting gene i 
nsh=zeros(n,1);%number of active shrna targeting gene i
ndes=zeros(n,1);%number of shrna designed to target gene i
v=zeros(n,1);%variance in mean log odds 
pt=pos_cnt>mad(pos_cnt(pos_cnt>0),1);%pt weeds out under sequenced reads in the treatment population
istar=zeros(n,1);%set to 1 if gene has any targeting shrna and zero otherwise
odr=odzp./odzn;%odds ratio
thr=1;%threshold for hairpin activity
po=pt&(odr>thr);%active hairpins
if size(uniq_locs,2)>size(uniq_locs,1),uniq_locs=uniq_locs';end
parfor i=1:n %for each gene compute the numbers of infected, sequenced, and active hairpins and the cumulative hairpin activity level
    uloc=uniq_locs==i;%hairpins targeting gene i
    ndes(i)=length(find(uloc));
    ntar(i)=length(find(uloc&pt));%sequenced hairpins targeting gene i
    pidx=find(uloc&po);%active hairpins targeting gene i
    nsh(i)=length(pidx);
    istar(i)=~isempty(pidx);
    if istar(i)
        %DerSimonian and Laird Random effects model
        yi=log((pc(pidx).*(nsum-nc(pidx)))./(nc(pidx).*(psum-pc(pidx))));
        wi=1./((pc(pidx)+nc(pidx))./(pc(pidx).*nc(pidx))+((nsum-nc(pidx))+(psum-pc(pidx)))./((nsum-nc(pidx)).*(psum-pc(pidx))));
        yw=sum(wi.*yi)/sum(wi);
        Q=sum(wi.*(yi-yw).^2);
        del=max(0,(Q-length(wi)-1)/(sum(wi)-sum(wi.^2)/sum(wi)));
        ws=1./(1./wi+del);
        mlodz(i)=sum(ws.*yi)/sum(ws);
        v(i)=1/sum(ws);
    else
        mlodz(i)=log(mean(odr(uloc&pt)));	
    end
end
mlodz(isnan(mlodz))=0;
%M1 gives the probability that gene i pareto dominates gene j in mean log
%odds, M2 likewise for number of targeting shRNA
M1=zeros(n,n);M2=zeros(n,n);
parfor i=1:n %for all genes compute all pairwise domination probabilities in both metrics 
    tmp1=zeros(1,n);
    tmp2=zeros(1,n);
    for j=1:i-1
        if (istar(i)|istar(j))&mlodz(i)>0&mlodz(j)>0
	    tmp1(j)=(mlodz(i)-mlodz(j))/(sqrt(v(i)+v(j)));
        %tmp2(j)=(nsh(i)-nsh(j))/(sqrt(nsh(i)+nsh(j)));
	    %d=ntar(i)/ntar(j)
	    tmp2(j)=sqrt(2)*(sqrt(nsh(i)+3/8)-sqrt(nsh(j)+3/8))/sqrt(sqrt(nsh(i)+3/8)+sqrt(nsh(j)+3/8));
        end
    end
    M1(i,:)=tmp1;M2(i,:)=tmp2;
end
M1=M1-tril(M1)';M2=M2-tril(M2)';
%diag(diag(A)) creates a diagonal matrix from the main diag of A
M1=.5*(1+tanh(M1/(.8*sqrt(2))));M1=M1-diag(diag(M1));
M2=.5*(1+tanh(M2/(.8*sqrt(2))));M2=M2-diag(diag(M2));
pd=sum(M1.*M2,2);%sum of probabilities of domination
pnd=sum((1-M1).*(1-M2),2);%sum of probabilities of not dominating
prnk=pd+.5*(size(M1,1)-pd-pnd);%sum of the prob of domination and .5 the prob of not being comparable
prnk(nsh==0)=0;%although genes without active shRNA are allowed to compete in the ranking if their 
%mlodz is >0 (they may have all undersequenced shRNA with none the less positive log odds),
%they are disqualified from earning a rank greater than a gene with active
%shRNA, rather they will be ranked by mean log odds

%rank genes first by mean log odds, then by domination probability
[~,trnk]=sortrows(mlodz,-1);
[~,tidx]=sortrows(prnk(trnk),-1);
trnk=trnk(tidx);rnk=zeros(n,1);
rnk(trnk)=1:n;


