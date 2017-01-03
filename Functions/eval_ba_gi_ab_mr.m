function eval_ba_gi_ab_mr(pdata,mdata,dbh_lbounds)

% =========================================================================
% This function evaluates the size structure of
% basal area _ba
% growth increment _gi
% abundance _ab
% mortality rate _mr
%
% If census data is available, then we compare
% =========================================================================

common_constants;

ndbhc          = length(dbh_lbounds); % Number of dbh classes

% Part I. Process Census data if it exists

if(pdata.exists)
    
    
    npy = numel(pdata.years);
    nha = pdata.nha;
    
    ipy_inc  = 1:npy-1;   % Indices for plot years for growth increments
    ipy_mass = 2:npy;     % Index for plot years for mass and stemcount
    
    
    stemd          = zeros(ndbhc,npy-1);     % Stem density
    incr           = zeros(ndbhc,npy-1);     % Mass increment [Mg/ha/yr]
    dincr          = zeros(ndbhc,npy-1);
    bincr          = zeros(ndbhc,npy-1);
    bmass          = zeros(ndbhc,npy-1);     % Mass [Mg]
    basal          = zeros(ndbhc,npy-1);
    kappa          = zeros(ndbhc,npy-1);     % Mortality counts
    
    % When calculating biomass, use all available data from this year, but
    % remove the data if there is weirdness in the growth increment, but
    % if there is no data from the previous year, don't calculate calculate
    % the delta
    
    
    ipc=0;
    for ip=ipy_mass
        
        ipc=ipc+1;
        deltayr  = pdata.years(ip)-pdata.years(ip-1);
        dbh      = pdata.dbh(:,ip);
        dbh_p    = pdata.dbh(:,ip-1);
        dddy     = (dbh-dbh_p)./deltayr;   % Change in D per change in year
        
        ids = find( (dbh>0 & dbh_p>0 & dddy>-0.5 & dddy<4.0) | ...
            (dbh>0 & dbh_p<0) );
        
        dbh   = dbh(ids);
        mass   = kgtomg*exp(-0.37+0.333.*log(dbh)+0.933.*log(dbh).^2.0-0.122.*log(dbh).^3.0);
        nst = numel(mass);
        
        tot_agb(ipc) = sum(mass)./nha;
        for idb=1:ndbhc-1
            dbh_lb = dbh_lbounds(idb);   % lower bound of bin
            dbh_ub = dbh_lbounds(idb+1); % upper bound of bin
            ids=find(dbh>dbh_lb & dbh<dbh_ub);
            stemd(idb,ipc)  = numel(ids)/nha;
            bmass(idb,ipc)  = sum(mass(ids))/nha;
            basal(idb,ipc)  = sum(pi*((dbh(ids)/100).^2)/4)/nha;
        end
        
        % DO THE LAST INDEX (UNBOUNDED)
        idb = ndbhc;
        dbh_lb = dbh_lbounds(idb);   % lower bound of bin
        ids=find(dbh>dbh_lb);
        stemd(idb,ipc)  = numel(ids)/nha;
        bmass(idb,ipc)  = sum(mass(ids))/nha;
        basal(idb,ipc)  = sum(pi*((dbh(ids)/100).^2)/4)/nha;
        
        
    end
    
    ipc=0;
    for ip=ipy_inc
        
        ipc=ipc+1;
        deltayr  = pdata.years(ip+1)-pdata.years(ip);
        dbh      = pdata.dbh(:,ip);
        dbh_n    = pdata.dbh(:,ip+1);                  %DBH of next year's data
        ba       = pi*((dbh/100).^2)/4;                % BA in m^2;
        ba_n     = pi*((dbh_n/100).^2)/4;              % next year's BA
        dddy     = (dbh_n-dbh)./deltayr;   % Change in D per change in year
        bady     = (ba_n-ba)./deltayr;
        
        % In calculating increments, use only data where both the current
        % and next year are positive, and are within reasonable limits
        % =====================================================================
        ids = find(dbh>0 & dbh_n>0  & dddy>-0.5 & dddy<4.0);
        
        % Ids of tree growing (as long as the second year is there),
        % this includes new recruits that were not there the first year.
        %	ids = find(dbh_n>0 & dddy>-0.5 & dddy<4.0 );
        %    ids_new = find(~(dbh(ids)>0.0001));
        %    dbh(ids(ids_new)) = 0;
        
        dbh   = dbh(ids);
        dbh_n = dbh_n(ids);
        dddy  = dddy(ids);
        bady  = bady(ids);
        
        mass   = kgtomg*exp(-0.37+0.333.*log(dbh)+0.933.*log(dbh).^2.0-0.122.*log(dbh).^3.0);
        mass_n = kgtomg*exp(-0.37+0.333.*log(dbh_n)+0.933.*log(dbh_n).^2.0-0.122.*log(dbh_n).^3.0);
        delb   = (mass_n-mass)/deltayr;
        
        tot_inc(ipc) = sum(delb)./nha;
        
        
        for idb=1:ndbhc-1
            
            dbh_lb = dbh_lbounds(idb);   % lower bound of bin
            dbh_ub = dbh_lbounds(idb+1); % upper bound of bin
            ids=find(dbh>dbh_lb & dbh<dbh_ub);
            incr(idb,ipc)   = sum(delb(ids))/nha;
            dincr(idb,ipc)  = mean(dddy(ids));
            bincr(idb,ipc)  = sum(bady(ids))/nha;
            
        end
        % Last un-bounded size class
        idb = ndbhc;
        dbh_lb = dbh_lbounds(idb);   % lower bound of bin
        ids=find(dbh>dbh_lb);
        incr(idb,ipc)   = sum(delb(ids))/nha;
        dincr(idb,ipc)  = mean(dddy(ids));
        bincr(idb,ipc)  = sum(bady(ids))/nha;
        
        
        % For the fluxes, we remove the constraints on positives and such
        % So re-read the data without the positive filters.  We do want to
        % remove crazy measurements, so keep the <4 >-0.5 filters.  But, we
        % want to keep those stems that were non-existant at the prior or
        % died.  So we want to keep it if the difference is >100 or <100
        % too...
        % =================================================================
        deltayr     = pdata.years(ip+1)-pdata.years(ip);
        dbh     = pdata.dbh(:,ip);
        dbh_n   = pdata.dbh(:,ip+1);
        dddy = (dbh_n-dbh)./deltayr;
        ids  = find( (dddy<4.0 | dddy>100) & (dddy>-0.5 | dddy<-100));
        
        dbh   = dbh(ids);
        dbh_n = dbh_n(ids);   %dbh_n is the dbh of the next census
                              % we use this to identify dead trees
        
        
        for idb=1:ndbhc
            
            dbh_lb = dbh_lbounds(idb);   % lower bound of bin
            if(idb<ndbhc)
                dbh_ub = dbh_lbounds(idb+1); % upper bound of bin
            else
                dbh_ub = Inf;
            end
            
            % Determine the ids of trees that died over that period
            n_mort                  = length(find(dbh>dbh_lb & dbh<dbh_ub & dbh_n<0));
            n_cnt                   = length(find(dbh>dbh_lb & dbh<dbh_ub));
            kappa(idb,ipc) = -log(1-(n_mort./n_cnt))/deltayr;
            
        end
        
    end
    
    
end

% =========================================================================
% Part II. Process Model Data
% =========================================================================

if(pdata.exists)
    
    immy = ed_census_overlap_indices(pdata.years,[mdata.date(:,1).year]);
    
else
    npy = numel(mdata.date(:,1).year);
    immy = 1:npy;
end

nmy     = numel(immy);

stemd_ed       = zeros(ndbhc,nmy);     % Stem density
incr_ed        = zeros(ndbhc,nmy);     % Mass increment [Mg/ha/yr]
dincr_ed       = zeros(ndbhc,nmy);
bincr_ed       = zeros(ndbhc,nmy);
bmass_ed       = zeros(ndbhc,nmy);     % Mass [Mg]
basal_ed       = zeros(ndbhc,nmy);
mrate_ed        = zeros(ndbhc,nmy);
mrate1_ed       = zeros(ndbhc,nmy);
mrate2_ed       = zeros(ndbhc,nmy);
mrate3_ed       = zeros(ndbhc,nmy);
mrate4_ed       = zeros(ndbhc,nmy);
mrate5_ed       = zeros(ndbhc,nmy);

bmrate_ed      = zeros(ndbhc,nmy);
tot_inc_ed     = zeros(nmy,1);
tot_agb_ed     = zeros(nmy,1);

otw = 1/12; %"over-twelve"
m2ha = 10000;  % sq meters to ha

iyc=0;
for iy  = immy
    
    iyc = iyc+1;
    year_ed(iyc) = mdata.date(iy,6).year;
    for im=1:12
        cted   = mdata.date(iy,im);
        npatch = cted.npatch;
        
        nco=0;
        for ipa = 1:npatch
            nco = nco+length(cted.pa(ipa).ed_codbh);
        end
        ed_raw_bmass  = zeros(nco,1);
        ed_raw_ba     = zeros(nco,1);
        ed_raw_delb   = zeros(nco,1);
        ed_raw_delba  = zeros(nco,1);
        ed_raw_deld   = zeros(nco,1);
        ed_raw_n      = zeros(nco,1);
        ed_raw_dbh    = zeros(nco,1);
        ed_raw_pft    = zeros(nco,1);
        ed_raw_mr     = zeros(nco,6);
        
        ic = 0;
        
        for ipa = 1:npatch
            afrac = cted.pa(ipa).afrac;
            nc = numel(cted.pa(ipa).ed_codbh);
            for ico=1:nc
                ic = ic+1;
                ed_raw_bmass(ic)  = kgtomg*C2B*cted.pa(ipa).ed_coagb(ico);
                ed_raw_dbh(ic)    = cted.pa(ipa).ed_codbh(ico);
                ed_raw_delb(ic)   = kgtomg*C2B*cted.pa(ipa).ed_dagb(ico);
                ed_raw_deld(ic)   = cted.pa(ipa).ed_ddbh(ico);
                ed_raw_n(ic)      = m2ha*afrac*cted.pa(ipa).ed_con(ico);
                ed_raw_ba(ic)     = pi*((cted.pa(ipa).ed_codbh(ico)/100.0).^2)/4;
                ed_raw_delba(ic)   = pi*(((cted.pa(ipa).ed_codbh(ico)+ed_raw_deld(ic))/100.0).^2)/4-ed_raw_ba(ic);
                ed_raw_pft(ic)    = cted.pa(ipa).ed_copft(ico);
                ed_raw_mr(ic,1)   = cted.pa(ipa).ed_mr(1,ico);
                ed_raw_mr(ic,2)   = cted.pa(ipa).ed_mr(2,ico);
                ed_raw_mr(ic,3)   = cted.pa(ipa).ed_mr(3,ico);
                ed_raw_mr(ic,4)   = cted.pa(ipa).ed_mr(4,ico);
                ed_raw_mr(ic,5)   = cted.pa(ipa).ed_mr(5,ico);
                ed_raw_mr(ic,6)   = sum([cted.pa(ipa).ed_mr(:,ico)]);
            end
        end
        
        tot_agb_ed(iyc) = tot_agb_ed(iyc)+otw*sum(ed_raw_bmass.*ed_raw_n);
        tot_inc_ed(iyc) = tot_inc_ed(iyc)+otw*sum(ed_raw_delb.*ed_raw_n);
        
        for idb=1:ndbhc
            dbh_lb = dbh_lbounds(idb);   % lower bound of bin
            
            if(idb<ndbhc)
                dbh_ub = dbh_lbounds(idb+1); % upper bound of bin
            else
                dbh_ub = Inf;
            end
            
            ids=find(ed_raw_dbh>dbh_lb & ed_raw_dbh<dbh_ub);
            
            stemd_ed(idb,iyc)  = stemd_ed(idb,iyc)+otw*sum(ed_raw_n(ids));
            incr_ed(idb,iyc)   = incr_ed(idb,iyc)+otw*sum(ed_raw_delb(ids).*ed_raw_n(ids));
            dincr_ed(idb,iyc)  = dincr_ed(idb,iyc)+otw*sum(ed_raw_deld(ids).*ed_raw_n(ids));
            bincr_ed(idb,iyc)  = bincr_ed(idb,iyc)+otw*sum(ed_raw_delba(ids).*ed_raw_n(ids));
            bmass_ed(idb,iyc)  = bmass_ed(idb,iyc)+otw*sum(ed_raw_bmass(ids).*ed_raw_n(ids));
            basal_ed(idb,iyc)  = basal_ed(idb,iyc)+otw*sum(ed_raw_ba(ids).*ed_raw_n(ids));  %m2/plant   * plant/m2  = m2/m2;
            mrate_ed(idb,iyc)   =  mrate_ed(idb,iyc)+otw*sum(ed_raw_n(ids).*ed_raw_mr(ids,6))./sum(ed_raw_n(ids));
            mrate1_ed(idb,iyc)  = mrate1_ed(idb,iyc)+otw*sum(ed_raw_n(ids).*ed_raw_mr(ids,1))./sum(ed_raw_n(ids));
            mrate2_ed(idb,iyc)  = mrate2_ed(idb,iyc)+otw*sum(ed_raw_n(ids).*ed_raw_mr(ids,2))./sum(ed_raw_n(ids));
            mrate3_ed(idb,iyc)  = mrate3_ed(idb,iyc)+otw*sum(ed_raw_n(ids).*ed_raw_mr(ids,3))./sum(ed_raw_n(ids));
            mrate4_ed(idb,iyc)  = mrate4_ed(idb,iyc)+otw*sum(ed_raw_n(ids).*ed_raw_mr(ids,4))./sum(ed_raw_n(ids));
            mrate5_ed(idb,iyc)  = mrate5_ed(idb,iyc)+otw*sum(ed_raw_n(ids).*ed_raw_mr(ids,5))./sum(ed_raw_n(ids));
            bmrate_ed(idb,iyc) = bmrate_ed(idb,iyc)+otw*sum(ed_raw_bmass(ids).*ed_raw_n(ids).*ed_raw_mr(ids,6));
            
        end
        
    end
    
    
end

% Normalize the rate on the diameter increment, it was weighted by stem
% density.
dincr_ed = dincr_ed./stemd_ed;
mean_mort_ed = mean(mrate_ed,2);   %
mean_mort_ed1 = mean(mrate1_ed,2); %
mean_mort_ed2 = mean(mrate2_ed,2); %
mean_mort_ed3 = mean(mrate3_ed,2); %
mean_mort_ed4 = mean(mrate4_ed,2); %
mean_mort_ed5 = mean(mrate5_ed,2); %


% Prepare the arrays for the figures depending on if census data is available

if(pdata.exists)
    
    stemd_both = [mean(stemd_ed,2)';mean(stemd,2)'];
    incr_both  = [mean(dincr_ed,2)';mean(dincr,2)'];
    basal_both = [mean(basal_ed,2)';mean(basal,2)'];
    mean_mort  = [mean(mrate_ed,2)';mean(kappa,2)'];
    
else
    
    stemd_both = mean(stemd_ed,2)';
    incr_both  = mean(dincr_ed,2)';
    basal_both = mean(basal_ed,2)';
    mean_mort  = mean(mrate_ed,2)';
    
end

set(0,'DefaultAxesFontSize',fasz);

[dbh_ticks,dbh_label,dbh_points] = dbh_axis_vectors(dbh_lbounds);


% BASAL AREA
% =========================================================================

figure;
set(gcf,'Units','Inches','Position',[1.0,1.0,5,4],'Color','w')
h1=plot(dbh_points,basal_both);
set(gca,'FontSize',fasz,'XTick',dbh_ticks);
set(h1(1),'Marker','o','MarkerFaceColor',cd.drk_rd,'MarkerEdgeColor','k','Color',cd.lgt_rd,'MarkerSize',8);
if(pdata.exists)
    set(h1(2),'Marker','o','MarkerFaceColor',cd.lgt_rd,'MarkerEdgeColor','k','Color',cd.lgt_rd,'MarkerSize',8);
end
ylabel('Basal Area [m^2 ha^{-1}]');
set(gca,'XLim',[min(dbh_ticks),max(dbh_ticks)]);
xlabel('DBH [cm]','FontSize',fasz);
box on;
grid on;
sumstr = sprintf('Totals \nED2: %3.1f',sum(basal_both(1,:)));
if(pdata.exists)
    sumstr = sprintf('%s \nField: %3.1f [m^2 ha^{-1}]',sumstr,sum(basal_both(2,:)));
    legend('ED2','Census','Location','NorthEast');
end
gtext(sumstr,'FontSize',fasz-2,'BackgroundColor','w');


% DIAMETER INCREMENT
% =========================================================================
figure;
set(gcf,'Units','Inches','Position',[1.5,1.5,5,4],'Color','w')

h1=plot(dbh_points,incr_both);
set(gca,'FontSize',fasz,'XTick',dbh_ticks);
set(h1(1),'Marker','o','MarkerFaceColor',cd.drk_rd,'MarkerEdgeColor','k','Color',cd.lgt_rd,'MarkerSize',8);
if(pdata.exists)
    set(h1(2),'Marker','o','MarkerFaceColor',cd.lgt_rd,'MarkerEdgeColor','k','Color',cd.lgt_rd,'MarkerSize',8);
    legend('ED2','Census','Location','NorthEast');
end
ylabel('Diameter Increment [cm yr^{-1}]');
xlabel('DBH [cm]','FontSize',fasz);
set(gca,'XLim',[min(dbh_ticks),max(dbh_ticks)]);
box on;
grid on;



% ABUNDANCE
% =========================================================================
figure;
set(gcf,'Units','Inches','Position',[2.0,2.0,5,4],'Color','w')
h1=semilogy(dbh_points,stemd_both);
set(gca,'FontSize',fasz,'XTick',dbh_ticks);
set(h1(1),'Marker','o','MarkerFaceColor',cd.drk_rd,'MarkerEdgeColor','k','Color',cd.lgt_rd,'MarkerSize',8);
if(pdata.exists)
    set(h1(2),'Marker','o','MarkerFaceColor',cd.lgt_rd,'MarkerEdgeColor','k','Color',cd.lgt_rd,'MarkerSize',8);
    legend('ED2','Census','Location','NorthEast');
end
ylabel('Abundance [ha^{-1}]');
xlabel('DBH [cm]','FontSize',fasz);
set(gca,'XLim',[min(dbh_ticks),max(dbh_ticks)]);
box on;
grid on;


% MORTALITY COMPARISON
% =========================================================================
figure;
set(gcf,'Units','Inches','Position',[2.5,2.5,5,4],'Color','w')

plth = plot(dbh_points,mean_mort);
set(plth(1),'Marker','o','MarkerFaceColor',cd.drk_rd,'MarkerEdgeColor','k','Color',cd.lgt_rd,'MarkerSize',8);
if(pdata.exists)
    set(plth(2),'Marker','o','MarkerFaceColor',cd.lgt_rd,'MarkerEdgeColor','k','Color',cd.lgt_rd,'MarkerSize',8);
    legend('ED2','Census','Location','NorthEast');
end
set(gca,'XLim',[min(dbh_ticks),max(dbh_ticks)]);
set(gca,'XTickLabel',dbh_label,'FontSize',fasz,'XTick',dbh_ticks);
set(gca,'YLim',[0,1.1*max([max(mean_mort),max(mean_mort_ed)])]);
xlabel('DBH Class [cm]','FontSize',13);
ylabel('Mortality Rate [year^{-1}]','FontSize',13);
%legend('Field','ED2');
title('rate = -log(1-(\DeltaN/N))/\Deltat;','Interpreter','tex');
grid on;
box on;



% MORTALITY PARTITIONS
% =========================================================================
figure;
set(gcf,'Units','Inches','Position',[3.0,3.0,5,4],'Color','w')

plth = plot(dbh_points,mean_mort_ed1,dbh_points,mean_mort_ed2,dbh_points,mean_mort_ed3+mean_mort_ed5,dbh_points,mean_mort_ed4);

set(plth(1),'Marker','o','MarkerEdgeColor',[0.0 0.3 0.0],'MarkerFaceColor',[0.6 0.8 0.6],'Color',[0.6 0.8 0.6],'MarkerSize',8);
set(plth(2),'Marker','o','MarkerEdgeColor',[0.0 0.0 0.3],'MarkerFaceColor',[0.6 0.6 0.8],'Color',[0.6 0.6 0.8],'MarkerSize',8);
set(plth(3),'Marker','o','MarkerEdgeColor',[0.3 0.0 0.0],'MarkerFaceColor',[0.8 0.6 0.6],'Color',[0.8 0.6 0.6],'MarkerSize',8);
set(plth(4),'Marker','o','MarkerEdgeColor',[0.3 0.3 0.0],'MarkerFaceColor',[0.8 0.8 0.6],'Color',[0.8 0.8 0.6],'MarkerSize',8);

ylim([0 0.04]);
set(gca,'XLim',[min(dbh_ticks),max(dbh_ticks)]);
set(gca,'XTickLabel',dbh_label,'FontSize',fasz,'XTick',dbh_ticks);
lhan=legend('Senescence','Carbon Balance','Background','Frost');
set(lhan,'FontSize',fasz);
%title('ED2 Mortality Rate Breakdown');

xlabel('DBH Class [cm]','FontSize',fasz);
ylabel('ED Mortality Parts [year^{-1}]','FontSize',fasz);

grid on;
box on;





