function oidx=pareto_plot(prank,nsh,mlodz)
%function pareto_plot(prank,nsh,mlodz)
%
%IN: prank is the pareto gene ranking vector
%    nsh(i) is the number of active shrna targeting gene i
%    mlodz(i) is the mean log odds ratio of shrna targeting gene i

prank=max(prank)-prank;
pct=[.1 .2 .3 .4 .5 .6 .7 .8 .9 .91 .92 .93 .94 .95 .96 .97 .98 .99 1];
q=quantile(prank,pct);
pidx=(mlodz>0);
nidx=(nsh>0);
oidx=pidx&nidx;
nsh=nsh(oidx);
mlodz=mlodz(oidx);
prank=prank(oidx);
grp=ones(size(mlodz));ngrps=1;
idx=find(prank<=q(1));
grp(idx)=ngrps;ngrps=ngrps+1;
for i=1:length(pct)-1
    idx=find(prank<=q(i+1)&prank>q(i));
    grp(idx)=ngrps;ngrps=ngrps+1;
end
s='';
for i=1:length(unique(grp))-1,s=[s 'o'];end
s=[s '*'];
cc=cool(length(unique(grp))-1);
size(cc)
length(unique(grp))
gscatter(nsh,mlodz,grp,[cc;[1 0 0]],s,10,'off')