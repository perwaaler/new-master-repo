function [lbnd,rbnd] = search_range_dmin(loc_extr, ind)
% calculates the search range for finding the minimum of d_min(t).
% * loc_extr contains approximations of minima and maxima of d_min(t)
% * ind is a vector of indeces for positions that have been predicted

% n is the number of predicted positions
n = length(ind);
loc_min = loc_extr.loc_min;
loc_max = loc_extr.loc_max;

if isempty(loc_min)
    i_min = 0;
else
    i_min = loc_min(1);
end

if length(loc_min)==2
    % find index of bounds
    l_nbhr = ind( max([i_min-1 ,1]) );
    r_nbhr = ind( ind==loc_max);
    
    lbnd = ind2time(l_nbhr);
    rbnd = ind2time(r_nbhr);
    
elseif length(loc_extr.loc_min)==1
    % find index of bounds
    l_nbhr = ind( max([i_min-1 ,1]) );
    r_nbhr = ind( min([i_min+1 ,n]) );
    
    lbnd = ind2time(l_nbhr);
    rbnd = ind2time(r_nbhr);
else
    lbnd = 0;
    rbnd = ind2time(n);
end

end