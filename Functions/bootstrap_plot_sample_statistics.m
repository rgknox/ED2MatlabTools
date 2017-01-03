
function bootstrap_plot_sample_statistics(pdata)

common_constants;

nha     = pdata.nha;
kgtomgh = kgtomg/nha;  % kgtomg is a common constant
ntrees  = size(pdata.dbh,1);
npy     = numel(pdata.years);


sampsz  = zeros(ntrees,1);
sampcnt = zeros(npy,1);
samperr = zeros(npy,1);

srtdbh  = NaN*zeros(npy,3);
srtbms  = NaN*zeros(npy,3);

setdbh  = zeros(npy,5000);
setbms  = zeros(npy,5000);

nbsamp=1000;
bssamp=zeros(nbsamp,1);

display('============================================');
display('Bootstrapping sample statistics on plot data');
display('This will take a while');


for it=1:ntrees
    dbhvec = pdata.dbh(it,:);
    
    % Determine the continuous vector of good data
    
    ids=find(dbhvec>-100);
    nids=numel(ids);
    
    if(nids>1)
        
        dbhvec = dbhvec(ids);
        bmsvec = kgtomgh*exp(-0.37+0.333*log(dbhvec)+0.933*log(dbhvec).^2.0-0.122*log(dbhvec).^3.0);
        
        % Check for gaps?
        if (sum(ids)==sum(ids(1):ids(end)))
            
            sampsz(it) = length(dbhvec);
            
                %                    display(sprintf('%d of %d',it,ntrees)); %#ok<UNRCH>
                if(sampsz(it)>1)
                    for is=1:nbsamp
                        % Generate a set of random id's
                        bsids =ceil(sampsz(it)*rand(sampsz(it),1));
                        bssamp(is) = mean(dbhvec(bsids));
                    end
                    % Determine the variance in the boot-strap means
                    samperr(nids) = samperr(nids)+sqrt(mean((bssamp-mean(bssamp)).^2));
                    sampcnt(nids) = sampcnt(nids)+1;
                    setdbh(nids,sampcnt(nids)) = mean(dbhvec);
                    setbms(nids,sampcnt(nids)) = mean(bmsvec);
                end
           
        end %if continuous
    end %if(nids>1)
end

for k=2:npy
    nk = sampcnt(k);
    if(nk>10)
        nk10 = round(0.10*nk);
        nk90 = round(0.90*nk);
        nk50 = round(0.5*nk);
        sortdbh = sort(setdbh(k,1:nk));
        sortbms = sort(setbms(k,1:nk));
        
        srtdbh(k,1) = sortdbh(nk10);
        srtdbh(k,2) = sortdbh(nk50);
        srtdbh(k,3) = sortdbh(nk90);
        
        srtbms(k,1) = sortbms(nk10);
        srtbms(k,2) = sortbms(nk50);
        srtbms(k,3) = sortbms(nk90);
    end
    
end

samperr=samperr./sampcnt;
set(0,'DefaultAxesFontSize',fasz);


figure(1);
set(gcf,'Units','Inches','Position',[2.0,2.0,7,4],'Color','w')
subplot(1,2,1);
h1=plot(1:npy,samperr);
set(h1(1),'Marker','o','MarkerFaceColor',cd.lgt_bu,'MarkerEdgeColor',cd.drk_bu,'Color',cd.lgt_bu,'MarkerSize',8);
xlabel('Measurements per Tree (Max is Number of Census)');
xlim([1,9]);
ylabel('\epsilon');
box on;
grid on;
% Plot the mean dbh associated with sample size
subplot(1,2,2);
for i=2:npy
    line([i,i],[srtdbh(i,1),srtdbh(i,3)],'Color',cd.lgt_bu);
end
hold on;
h2=plot(1:npy,srtdbh);
set(h2(2),'Marker','o','MarkerFaceColor',cd.lgt_bu,'MarkerEdgeColor',cd.drk_bu,'Color',cd.lgt_bu,'MarkerSize',8);
set(h2(1),'Marker','o','MarkerFaceColor',cd.lgt_gn,'MarkerEdgeColor',cd.drk_gn,'MarkerSize',8,'LineStyle','none');
set(h2(3),'Marker','o','MarkerFaceColor',cd.lgt_pu,'MarkerEdgeColor',cd.drk_pu,'MarkerSize',8,'LineStyle','none');

% Make some lines


ylim([0,40]);

xlabel('Measurements Per Tree');
xlim([1,9]);
ylabel('Diameter Percentiles');
box on;
grid on;

gtext('10th','FontSize',12);
gtext('50th','FontSize',12);
gtext('90th','FontSize',12);


