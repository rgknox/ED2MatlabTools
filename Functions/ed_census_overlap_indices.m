
function [ids] = ed_census_overlap_indices(pdata_years,mdata_years)

   
    % Perform a matching on time indices, if there is census data available
    % then use it to constrain the time points, otherwise just use all of the
    % model data.
    
    imm1=find(mdata_years>=pdata_years(1),1,'first');
    imm2=find(mdata_years<=pdata_years(end),1,'last');

    ids= imm1:imm2;  % The indices of the relevant years