function plot_fluxes(dset,varnames,varunits,phrs,fluxes_img,visible)

global fasz;

ntmp = length(dset(1).flux_tp(1,:));
nfigs = ceil(ntmp/4);

fasz_l=fasz+1;
ip=0;

for ifig=1:nfigs
    
    figure('visible',visible);
    set(gcf,'PaperPositionMode','manual',...
        'Units','inches','Position',[0.25+0.1*ifig 0.25+0.1*ifig 8.5 6]);
    
    bx = 0.09;mx = 0.07;  % Set the bottom and middle buffer spaces
    by = 0.09;my = 0.07;  % Set the same spaces for the vertical direction
    
    dx = ((0.99-bx)./4)-mx;
    dy = ((0.94-by)./2.0)-my;
    
    ndsets = numel(dset);
    
    
    symbols = {'o','+','*','^'};
    sycolor = {'b','r','g','k'};
    
    for pos=0:3
        
        ip=ip+1;
        ax2 = axes;
        set(ax2,'Position',[bx+pos*dx+pos*mx by+my+dy dx dy],'FontSize',fasz_l);
        hold on;
        minxy = 1e10;
        maxxy = -1e10;
        for kd=1:ndsets
            plot(phrs,dset(kd).flux_td(:,ip),symbols{kd},'Color',sycolor{kd},'MarkerSize',8);
            minxy = min([minxy,min(dset(kd).flux_td(:,ip))]);
            maxxy = max([maxxy,max(dset(kd).flux_td(:,ip))]);
        end
        minxy=minxy-0.03*abs(maxxy);
        maxxy=maxxy+0.03*abs(maxxy);
        
        ylim([minxy maxxy]);
        xlim([1 24]);
        set(gca,'XTick',[6 12 18]);
        hold off;
        grid on;
        box on;
        if(ip==ntmp);ylabel('Hourly Means','FontSize',fasz_l);end;
        title(sprintf('%s\n%s',varnames{ip},varunits{ip}),'FontSize',fasz_l);
        ax3 = axes;
        set(ax3,'Position',[bx+pos*dx+pos*mx by dx dy],...
            'FontSize',fasz_l);
        hold on;
        minxy = 1e10;
        maxxy = -1e10;
        for kd=1:ndsets
            plot(1:12,dset(kd).flux_tm(:,ip),symbols{kd},'Color',sycolor{kd},'MarkerSize',8);
            minxy = min([minxy,min(dset(kd).flux_tm(:,ip))]);
            maxxy = max([maxxy,max(dset(kd).flux_tm(:,ip))]);
        end
        hold off;
        xlim([1 12]);
        set(gca,'XTick',[3 6 9],'XTicklabel',{'Mar','Jun','Sep'},'FontSize',fasz_l);
        grid on;
        box on;
        if(ip==ntmp);ylabel('Monthly Means','FontSize',fasz_l);end;
        
        
        ax4 = axes;
        set(ax4,'Position',[0.0 0.0 1.0 by],...
            'FontSize',fasz_l);
        axis off;
        
        lablist = '';
        for kd=1:ndsets
            lablist = strcat(lablist,sprintf(' %s=%s ',dset(kd).lab,symbols{kd}));
        end
        text(0.45,0.3,lablist);
        
        if(ip==ntmp) 
            break
        end
    end
    
    oldscreenunits = get(gcf,'Units');
    oldpaperunits = get(gcf,'PaperUnits');
    oldpaperpos = get(gcf,'PaperPosition');
    set(gcf,'Units','pixels');
    scrpos = get(gcf,'Position');
    newpos = scrpos/100;
    set(gcf,'PaperUnits','inches',...
        'PaperPosition',newpos)
    print('-dpng',fluxes_img,'-r200');
    drawnow
    set(gcf,'Units',oldscreenunits,...
        'PaperUnits',oldpaperunits,...
        'PaperPosition',oldpaperpos)
    
end

%print('-dpng','-r300',fluxes_img);
%print('-djpeg','-r300','flux.jpg');

