
function [grfrac] = fraction_of_growing_trees(pdata,mdata,dbh_lbounds)

% =========================================================================
% This routine calculates the fraction of trees that are growing
% The first calulations are stored in an array, indexed by
%
% [ size class, time interval ]
%
% Note that comparing the model output and census data is a little tricky
% Census data comes at irregular intervals, typically 1 year or greater
% while model output from ED2 in the form of "Q" files comes every month.
% Note also that growth fractions from census data is a comparison of
% the same tree from different points in time, while ED2 simply gives us
% a rate associated with the cohort, the rate coming from model dynamics.
%
% Bearing all this in mind, the time dimension in the grfrac[] arrays for
% model and for census have different meaning.
% =========================================================================

common_constants;

near_zero = 0.99*0.1;      % This is the threshold for measurement precision
                           % which is currently 1 mm.  For instance
                           % the calculation of growth fractions will
                           % determine if any growth is happening beyond
                           % this precision.

                           
ndbhc     = length(dbh_lbounds); % Number of dbh classes (centers)


otw = 1/12; %"over-twelve"


% Growth Fractions From Census Data
% =========================================================================


if(pdata.exists)
    
    npy = length(pdata.years);
    imaxyears = round(pdata.years(end)-pdata.years(1));
    
    grfrac         = zeros(ndbhc,imaxyears);
    grsamp         = zeros(ndbhc,imaxyears);

    % Mortality rates of stagnant and growing trees 
    grmort         = zeros(ndbhc,imaxyears-1);
    grmort_c       = zeros(ndbhc,imaxyears-1);
    stmort         = zeros(ndbhc,imaxyears-1);
    stmort_c       = zeros(ndbhc,imaxyears-1);
    
    ntrees  = size(pdata.dbh,1);
    ndbhmax = length(pdata.dbh(1,:));
    
    for it=1:ntrees
        
        dbhvec = pdata.dbh(it,:);
        
        % Determine the continuous vector of good data
        
        ids=find(dbhvec>-100);
        nids=numel(ids);
        
        if(nids>1)
            
            dbhvec = dbhvec(ids);
            yearvec = pdata.years(ids);
            
            % Check for gaps?
            if (sum(ids)==sum(ids(1):ids(end)))
                
                % Calculate the growth fractions over different integration windows
                % Index with idy idy=1 means there is one span, ie two dbh records
                % ======================================================================
                for idy = 1:npy-1
                    % Calculate how many spans we can do
                    nspan = nids-idy;
                    if nspan>0
                        for in=1:nspan
                            
                            idyr   = round(yearvec(in+idy)-yearvec(idy));
                            ddbhdy = (dbhvec(in+idy)-dbhvec(idy))./(yearvec(in+idy)-yearvec(idy));
                            mdbh   = mean([dbhvec(idy):dbhvec(in+idy)]);
                            idb    = find(dbh_lbounds<mdbh,1,'last');
                            
                            % How to tell if the next year is a dead year?
                            % Global index in pdata.dbh: gid
                            gid = (ids(1)+in+idy);
                            ltry=false;
                            if(gid<=ndbhmax)
                                ltry=true;
                                if(pdata.dbh(it,gid)<-50)
                                    lmort=true;
                                    
                                else
                                    lmort=false;
                                end
                            end
                            
                            % This tree is actively growing
                            if(ddbhdy>=near_zero) %#ok<BDSCI>
                                grfrac(idb,idyr)=grfrac(idb,idyr)+1;
                                
                                if(ltry)
                                    if(lmort)
                                        grmort(idb,idyr)=grmort(idb,idyr)+1;
                                    end
                                    grmort_c(idb,idyr)=grmort_c(idb,idyr)+1;
                                end
                                % This tree is stagnant
                            else
                                
                                if(ltry)
                                    if(lmort)
                                        stmort(idb,idyr)=stmort(idb,idyr)+1;
                                    end
                                    stmort_c(idb,idyr)=stmort_c(idb,idyr)+1;
                                end
                                
                                
                                
                            end
                            grsamp(idb,idyr)=grsamp(idb,idyr)+1;
                        end
                    end
                end
                
            end %if continuous
        end %if(nids>1)
    end
    
    grmort = grmort./grmort_c;
    stmort = stmort./stmort_c;
    grfrac = grfrac./grsamp;
    
end

% =========================================================================
%
% Estimate the fraction of growing plants in ED2
% Important Note: It can be argued that since ED2 does not track individual
% trees, that it is assumed that there is variability in growth rates
% within each cohort and that the cohort is simply showing the mean
% behavior of that plant group.  So in this sense, it can be argued that it
% is not a fair comparison to look at the fraction of plants growing
% against cohort data.  As cohorts will show much less variability, and
% will not capture non-growers quite as well.
% But anyway ... if you don't
% care about that kind of thing, the following will calculate a growth
% fraction assuming all plants in the cohort have the same growth rates.
%
% =========================================================================

if(pdata.exists)
       
    immy = ed_census_overlap_indices(pdata.years,[mdata.date(:,1).year]);
    
else
    npy = numel(mdata.date(:,1).year);
    immy = 1:npy;
end

grfrac_ed = zeros(ndbhc,numel(immy));     % Growing Fraction
ncount_ed = zeros(ndbhc,numel(immy));

% Loop through all of the time-points, calculate growth between the current
% and the previous time-point.  See if that growth, (in biomass), is larger
% than the projected biomass growth of measurement precision.

iyc=0;
for iy  = immy
    
    iyc = iyc+1;
    for im=1:12
        % Model data is saved in structures by time-stamp
        % mdata.dat(:,:).stuff
        % cted is a pointer to the current Q file's date
        
        cted   = mdata.date(iy,im);    % <= current structure
        
        npatch = cted.npatch;
        nco=0;
        for ipa = 1:npatch
            nco = nco+length(cted.pa(ipa).ed_codbh);
        end
        ed_raw_delb   = zeros(nco,1);
        ed_raw_n      = zeros(nco,1);
        ed_raw_dbh    = zeros(nco,1);
        
        ic = 0;
        for ipa = 1:npatch
            afrac = cted.pa(ipa).afrac;
            nc = numel(cted.pa(ipa).ed_codbh);
            for ico=1:nc
                ic = ic+1;
                ed_raw_dbh(ic)    = cted.pa(ipa).ed_codbh(ico);
                ed_raw_delb(ic)   = kgtomg*C2B*cted.pa(ipa).ed_dagb(ico);
                ed_raw_n(ic)      = m2ha*afrac*cted.pa(ipa).ed_con(ico);  % Number of trees per hectare
            end
        end
        
        for idb=1:ndbhc
            dbh_lb = dbh_lbounds(idb);   % lower bound of bin
            
            if(idb<ndbhc)
                dbh_ub = dbh_lbounds(idb+1); % upper bound of bin
            else
                dbh_ub = Inf;
            end
            
            ids=find(ed_raw_dbh>dbh_lb & ed_raw_dbh<dbh_ub);
            
            % Calculate the equivalent delb
            
            equiv_b1 = kgtomg*exp(-0.37+0.333*log(ed_raw_dbh(ids))+...
                0.933*log(ed_raw_dbh(ids)).^2.0-0.122*log(ed_raw_dbh(ids)).^3.0);
            
            equiv_b2 = kgtomg*exp(-0.37+0.333*log(ed_raw_dbh(ids)+near_zero) + ...
                0.933*log(ed_raw_dbh(ids)+near_zero).^2.0-0.122*log(ed_raw_dbh(ids)+near_zero).^3.0);
            
            near_zero_delb = equiv_b2-equiv_b1;
            
            ids_g=[];
            for i=1:numel(ids)
                if(ed_raw_delb(ids(i))>near_zero_delb(i))
                    ids_g=[ids_g,ids(i)];
                end
            end
            
            if(numel(ids_g)>0)
                grfrac_ed(idb,iyc) = grfrac_ed(idb,iyc)+otw*sum(ed_raw_n(ids_g))./sum(ed_raw_n(ids));
                ncount_ed(idb,iyc) = ncount_ed(idb,iyc) + 1;
            end            
        end
        
    end
    
    
end


% Make some plots
% Note that there may be years or diameter classes that do not have
% data.  We must filter these out, and also determine if we have just model
% data, or both model and plot (census) data.
% =========================================================================

if(pdata.exists)
    
    [~,nuyrs] = size(grfrac);   % (size,year)
    
 
        % Note that the census growth fractions are from a number of
        % different time intervals!!!!!
        % In this step we must identify the time intervals that are present
        % because grfrac() is a sparse matrix (probably)
        % idyr is the index of year intevals (delta year = dyr)
        usedyr = [];
        for idyr = 1:nuyrs
            if (any(grsamp(:,idyr)))       % grsamp(idb,idyr)
                usedyr = [usedyr,idyr]; %#ok<AGROW>
            end
        end
        
        % OVERRIDE FOR PICKING SPECIFIC INTER-CENSUS INTERVALS
        % Trim usedyr to the 1st and last intervals
        %
        usedyr = [usedyr(1),usedyr(2),usedyr(end)];
        
        nintervals  = numel(usedyr);
        grfrac_mean = NaN*zeros(nintervals+1,ndbhc);

        legend_str = {'ED2'};

        for idb = 1:ndbhc
            ivalid = find(ncount_ed(idb,:)>0);
            if(numel(ivalid)>0)
                grfrac_mean(1,idb) = mean(grfrac_ed(idb,ivalid),2);
            end
        end    
        
        for iusedyr=1:numel(usedyr)
            legend_str{iusedyr+1} = sprintf('Census \\Delta %i yr',usedyr(iusedyr));
            for idb = 1:ndbhc
                if(grsamp(idb,iusedyr+1)>0)
                    grfrac_mean(iusedyr+1,idb) = grfrac(idb,usedyr(iusedyr));
                end
            end
        end
        
else
    
       grfrac_mean = NaN*zeros(1,ndbhc);
       for idb = 1:ndbhc
            ivalid = find(ncount_ed(idb,:)>0);
            if(numel(ivalid)>0)
                grfrac_mean(1,idb) = mean(grfrac_ed(idb,ivalid),2);
            end
       end    
       legend_str = 'ED2';
end



set(0,'DefaultAxesFontSize',fasz);

if(length(grfrac_mean(:,1))>1)
    % The lower most gray is too white, add an index
    cmap = cbrewer('seq','Greys',length(grfrac_mean(:,1))+1);
else
    cmap = [0.5,0.5,0.5];
end

[dbh_ticks,dbh_label,dbh_points] = dbh_axis_vectors(dbh_lbounds);

figure;
p_han = plot(dbh_points,grfrac_mean);

for i=1:length(grfrac_mean(:,1))
    set(p_han(i),'Color',cmap(i+1,:),'Marker','o','MarkerSize',marksize,'MarkerFaceColor',cmap(i+1,:),'MarkerEdgeColor','k')
end

ylim([0,1]);
set(gca,'XTickLabel',dbh_label,'FontSize',13,'XTick',dbh_ticks);
xlabel('DBH [cm]','FontSize',13);
ylabel('Growth Fraction','FontSize',13)
box on;grid on;
l_han=legend(legend_str);
set(l_han,'Interpreter','tex','Location','SouthEast')
title(sprintf('Fraction of Plants Showing \n Growth Increment > %2.1f cm',near_zero));

display('Press The AnyKey');
pause;

