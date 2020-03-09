function out = coll_type(data_type)
% function that gets the collision type of that the data type estimates.
% 1=DAFEA 2=X 3=ttc_FEA 4=ttc_min ttc_min_EA=5 danger_FEA=6 dist_min_EA=7 dist_min=8

if data_type == 1 || data_type == 3 || data_type == 6
    out = 1;
elseif data_type == 5 || data_type == 7
    out = 2;
else
    out = 3;
end
end