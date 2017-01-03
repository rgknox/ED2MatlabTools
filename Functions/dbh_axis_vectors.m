function [dbh_ticks,dbh_label,dbh_points] = dbh_axis_vectors(dbh_lbounds)


% =========================================================================
%
% This simple function creates the vectors of strings and real numbers that
% help with plotting a dbh axis.  The lower dbh bins are provided by the
% user, this function does two things.  It creates an extra end-point entry
% for the purposes of spacing the last bin, it also creates stings for
% labeling the edges of the bins and gives a "+" to the infinite upper edge
% of the last bin.  It finally returns a point for the bin center where the
%
% =========================================================================

ndbh = length(dbh_lbounds);

for i=1:ndbh
    dbh_label{i} = sprintf('%d',dbh_lbounds(i)); %#ok<*AGROW>
    dbh_ticks(i) = dbh_lbounds(i);
end
dbh_label{i+1} = '+';
dbh_ticks(i+1) = dbh_ticks(i) + (dbh_lbounds(i)-dbh_lbounds(i-1));

for i=1:ndbh
    dbh_points(i) = mean(dbh_ticks(i:i+1));
end