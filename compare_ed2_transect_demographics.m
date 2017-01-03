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
%
%==========================================================================

clear all;
close all;

addpath('Functions');
addpath('/home/rgknox/local/MATLAB/cbrewer');

common_constants;   % Initialize various constants to memory

% Analysis Control Parameters
% =========================================================================

have_pdata = false;  % Do you have plot-data?

dbh_lbounds = [10,20,30,40,50,60];


% Provide some information about the transect data if you have it.
% This assumes that the transect data is in a CSV file.
% IMPORTANT: The last year in the CSV file will not be used. The assumption
% is that the last year does not have the benefit of growth filtering. This
% was suggested by Jeff Chambers. If you really want that last year for
% static quantities like number density, mass and basal area, then... I
% don't know, I guess you will have to muck with the code.
% =========================================================================

if(have_pdata)
    tranfiles = {'/home/rgknox/Sync/Projects/LDRD_MODEX_CO2/ManausCensus/Transecto_2013_v1_rkfilled_LO.csv', ...
        '/home/rgknox/Sync/Projects/LDRD_MODEX_CO2/ManausCensus/Transecto_2013_v1_rkfilled_NS.csv'};
    theaders = {'UA','N','Espécie','DAP 96','DAP 00','DAP 02','DAP 04','DAP 06','DAP 08','DAP 10','DAP 11','DAP 13'};
    transect_years = [1996,2000,2002,2004,2006,2008,2010,2011,2013];
    dbh_columns    = [4,   5,   6,   7,   8,   9,   10,  11,  12];
    ignore_last_year = true;
    csvform = '%d %f %s %f %f %f %f %f %f %f %f %f';
    nhectares = [5.0,5.0];  % The total number of hectares in each part
    mort_flag = -666;           % Flag indicating mortality
    dnexist_flag = -444;        % Flag 
    pdata = load_pdata_csv(tranfiles,csvform,transect_years,dbh_columns,nhectares,ignore_last_year);
    pdata.exists = true;
else
    pdata.exists = false;
end

% =========================================================================
% Load ED Model output, using "Q"-files
% =========================================================================


edqpfx = 'TestData/ZF2_TF.9PCNT/zf2_oxi_TF.9PCNT-Q-';


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
if(do_eval_ba_gi_ab_mr)
    eval_ba_gi_ab_mr(pdata,mdata,dbh_lbounds);
end







%if(false)
%    objective_hists_edplots(nf,tf,ted,xbine_dap,nyear)
%end


% Objective 2
% Compare the sorted integrations of ED2 and census biomass



if(false)
%%    objective_comped_sortbio(pdata,mdata)
end

return;


