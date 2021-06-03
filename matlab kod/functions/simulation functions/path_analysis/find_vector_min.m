function out = find_vector_min(vec)
% calculates the index of the local minima and maxima.

n = length(vec);
local_min = [];
local_max = [];

for i=1:n
    
    if     i==1
        ldiff = -1;
        rdiff = vec(i) - vec(i+1);
        
    elseif i==n
        ldiff = vec(i) - vec(i-1);
        rdiff = -1;
        
    else
        ldiff = vec(i) - vec(i-1);
        rdiff = vec(i) - vec(i+1);
    end
    
    if max([ldiff,rdiff])<0
        local_min(end+1) = i; %#ok<*AGROW>
    end
    if min([ldiff,rdiff])>0
        local_max(end+1) = i;
    end
end

out.loc_min = local_min;
out.loc_max = local_max;

end