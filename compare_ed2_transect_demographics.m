%==========================================================================
%
% compare_ed2_transect_demographics.m
%
% This script evaluates ED2 demographics output. If transect data is
% present, it will compare with the transect data.
%
% RGK 09-2013
% Updated 01-2014
% Updated 12-2016
% Updated 01-2017
%
%==========================================================================

clear all;
close all;

addpath('Functions');
addpath('/home/rgknox/local/MATLAB/cbrewer');

common_constants;   % Initialize various constants to memory

% Analysis Control Parameters
% =========================================================================

have_pdata = true;  % Do you have plot-data?

dbh_lbounds = [10,20,30,40,50];

ci_level = 0.95;

% Provide some information about the transect data if you have it.
% This assumes that the transect data is in a CSV file.
% IMPORTANT: The last year in the CSV file will not be used. The assumption
% is that the last year does not have the benefit of growth filtering. This
% was suggested by Jeff Chambers. If you really want that last year for
% static quantities like number density, mass and basal area, then... I
% don't know, I guess you will have to muck with the code.
% =========================================================================

if(have_pdata)
%    tranfiles = {'/home/rgknox/Sync/Projects/LDRD_MODEX_CO2/ManausCensus/Transecto_2013_v1_rkfilled_LO.csv', ...
%        '/home/rgknox/Sync/Projects/LDRD_MODEX_CO2/ManausCensus/Transecto_2013_v1_rkfilled_NS.csv'};
%    theaders = {'UA','N','Espécie','DAP 96','DAP 00','DAP 02','DAP 04','DAP 06','DAP 08','DAP 10','DAP 11','DAP 13'};
%    transect_years = [1996,2000,2002,2004,2006,2008,2010,2011,2013];
%    site_name = 'Manaus ZF2';
%    dbh_columns    = [4,   5,   6,   7,   8,   9,   10,  11,  12];
%    ignore_last_year = true;
%    csvform = '%d %f %s %f %f %f %f %f %f %f %f %f';
%    nhectares = [5.0,5.0];  % The total number of hectares in each part
%    mort_flag = -666;           % Flag indicating mortality
%    dnexist_flag = -444;        % Flag 
%    pdata = load_pdata_csv(tranfiles,csvform,transect_years,dbh_columns,nhectares,ignore_last_year);
%    pdata.exists = true;
%    pdata.site_name = site_name;

    
%    tranfiles = {'TestData/gyf_summary_2004_2013.csv'};
%    site_name = 'GYF';
%    theaders = {'Trans','Tag','Name','DBH2004','DBH2006','DBH2008','DBH2010','DBH2013'};
%    transect_years = [2004,2006,2008,2010,2013];
%    dbh_columns    = [4,   5,   6,   7,   8];
%    ignore_last_year = true;
%    csvform = '%d %d %s %f %f %f %f %f';
%    nhectares = [4.91];  % The total number of hectares in each part
%    mort_flag = -666;           % Flag indicating mortality
%    dnexist_flag = -666;        % Flag 
%    pdata = load_pdata_csv(tranfiles,csvform,transect_years,dbh_columns,nhectares,ignore_last_year);
%    pdata.exists = true;
%    pdata.site_name = site_name;


    tranfiles = {'TestData/s67_summary_1999_2011.csv'};
    site_name = 'TNF';
    theaders = {'Trans','Tag','Name','DBH1999','DBH2001','DBH2005','DBH2008','DBH2009','DBH2010','DBH2011'};
    transect_years = [1999,2001,2005,2008,2009,2010,2011];
    dbh_columns    = [4,   5,   6,   7,   8,   9,   10];
    ignore_last_year = true;
    csvform = '%s %d %s %f %f %f %f %f %f %f';
    nhectares = [3.12];  % The total number of hectares in each part
    mort_flag = -666;           % Flag indicating mortality
    dnexist_flag = -666;        % Flag 
    pdata = load_pdata_csv(tranfiles,csvform,transect_years,dbh_columns,nhectares,ignore_last_year);
    pdata.exists = true;
    pdata.site_name = site_name;
    
    
%    trans,tag,name,dbh.1999,dbh.2001,dbh.2005,dbh.2008,dbh.2009,dbh.2010,dbh.2011
    
else
    pdata.exists = false;
end

% =========================================================================
% Load ED Model output, using "Q"-files
% =========================================================================


%edqpfx = 'TestData/ZF2_TF.9PCNT/zf2_oxi_TF.9PCNT-Q-';
%edqpfx = 'TestData/guyaflux/tgyf_yra-10_iphen-01_stext06_tfall140-Q-';
edqpfx = 'TestData/tapajos/ts67_yra-10_iphen-01_stext16_tfall140-Q-';

[mdata]=load_ed_mdata(edqpfx);


% Objective 1
% log(M) log(G) scatter plots like in Stephenson
%if(true)
%    objective_scatter_edsite(pdata,mdata);
%end

% =========================================================================
% Growth Fractions
% =========================================================================

do_growth_fractions = false;
if(do_growth_fractions)
    [grfrac]=fraction_of_growing_trees(pdata,mdata,dbh_lbounds);
end

% Sampling Statistics (this is very experimental, probably not worth your
% (time, as even I can't understand what I was trying to do :0 )
% =========================================================================

do_bootstrap_plot_sample_statistics = false;
if(do_bootstrap_plot_sample_statistics)
    bootstrap_plot_sample_statistics(pdata);
end

% =========================================================================
% Mean estimated basala area, growth increment, abundance and mortality
% rate, discretized by size class.
% =========================================================================

do_eval_ba_gi_ab_mr = true;
do_cis = true;
if(do_eval_ba_gi_ab_mr)
    eval_ba_gi_ab_mr(pdata,mdata,dbh_lbounds,do_cis,ci_level);
end

%if(false)
%    objective_hists_edplots(nf,tf,ted,xbine_dap,nyear)
%end


% Objective 2
% Compare the sorted integrations of ED2 and census biomass



%if(false)
%    objective_comped_sortbio(pdata,mdata)
%end

return;


