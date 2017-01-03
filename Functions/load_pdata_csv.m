function pdata = load_pdata_csv(tranfiles,csvform,transect_years,dbh_columns,nhectares,ignore_last_year)

% =========================================================================
% pdata, the main structure returned by this function: "plot-data"
% =========================================================================

if(size(dbh_columns)~=size(transect_years))
    display('The vector of column indices is not the same size as their');
    display('Associated year');
    pause;
    return;
end


cols_dbh         = dbh_columns;
nf               = length(tranfiles);
pdata.nha        = sum(nhectares);
if(ignore_last_year)
    nsyr             = length(dbh_columns)-1; % DONT INCLUDE THE LAST YEAR
else
    nsyr             = length(dbh_columns); % DO INCLUDE THE LAST YEAR
end

pdata.nyears     = nsyr;
pdata.years      = transect_years(1:nsyr);

% Open the plot data files, and first see how many rows we are dealing with
% The assumption is that the different files represent different parts of
% the plot of interest, such as different transects.
% =========================================================================

nrow=0;
for i=1:nf
    fid=fopen(tranfiles{i});
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        nrow=nrow+1;
    end
    fclose(fid);
end


% =========================================================================
% Initialize and load from a CSV file
% =========================================================================

pdata.tran  = zeros(nrow,1);            % Points back to the transect
pdata.ua    = zeros(nrow,1);            % This is the cluster number?
pdata.N     = zeros(nrow,1);            % Not sure
pdata.dbh   = zeros(nrow,nsyr);         % Diameter at breast height
pdata.dbh_n = zeros(nrow,nsyr);         % DBH of next year
pdata.bmass = zeros(nrow,nsyr);         % Allometry estimated biomass (Chambers)
pdata.delb  = zeros(nrow,nsyr);         % Growth increment



% Check data size and compare with the header names
irow = 0;  % Counter for all lines

for i=1:nf
    fid=fopen(tranfiles{i});
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        irow=irow+1;
        cline=textscan(tline,csvform,'delimiter',',');
        
        
        pdata.tran(irow)    = i;
        pdata.ua(irow)      = cline{1}(1);
        pdata.N(irow)       = cline{2}(1);
        for iy=1:nsyr
            
            % Do we have positive biomass (ie not dead or weird)
            % Then calculate a growth incrememnt
            
            dbh                = cline{cols_dbh(iy)}(1);
            pdata.dbh(irow,iy) = dbh;

        end
        
    end
    fclose(fid);
end
clear cols_dbh;

display('Finished loading census data');
