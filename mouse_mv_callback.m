function mouse_mv_callback(ax,pdata)
   %get current point and find closest data point in scatter plot
   if isempty(ax),ax=gca;end
   pt=get(ax,'currentpoint');x=round(pt(1,1));y=pt(1,2);
   xidx=find(x==pdata.xdat);
   if length(xidx)==0,return;end
   yd = (y-pdata.ydat);
   [~,yidx] = min(abs(yd(xidx)));
   idx=xidx(yidx);
   %set the current gene view values from the screen data
   try
   h=pdata.handles;
   set(h.gene_symbol,'String',pdata.lbls{idx});
   set(h.rank,'String',num2str(pdata.prank(idx)));
   if pdata.fdr(idx)<=0.05,
       set(h.screen_fdr_textbox,'String',num2str(pdata.fdr(idx)),'ForegroundColor','r');
   else
       set(h.screen_fdr_textbox,'String',num2str(pdata.fdr(idx)),'ForegroundColor','k');
   end
   set(h.nsh_textbox,'String',num2str(pdata.xdat(idx)));
   set(h.ntar_textbox,'String',num2str(pdata.ntar(idx)));
   if isempty(pdata.tt)|pdata.tt(idx)==-1
        set(h.ttest_textbox,'String','NA');
        set(h.fc_textbox,'String','NA');    
   else
       if pdata.tt(idx)<=0.05
        set(h.ttest_textbox,'String',num2str(pdata.tt(idx)),'ForegroundColor','r');
       else
        set(h.ttest_textbox,'String',num2str(pdata.tt(idx)),'ForegroundColor','k');
       end
        set(h.fc_textbox,'String',num2str(pdata.fc(idx)));
    end
   set(h.mlodz_textbox,'String',num2str(pdata.ydat(idx)));
   set(h.net_cent_textbox,'String',num2str(pdata.net_cent_comb(idx)));
   set(h.genetic_textbox,'String',num2str(pdata.net_cent_genetic(idx)));
   set(h.protein_textbox,'String',num2str(pdata.net_cent_protein(idx)));
   set(h.coloc_textbox,'String',num2str(pdata.net_cent_coloc(idx)));
   set(h.pathway_textbox,'String',num2str(pdata.net_cent_pathway(idx)));
   catch me
       return;
   end
end