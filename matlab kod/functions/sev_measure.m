function out = sev_measure(data_type)
% translates data type to severity measure

if data_type == 1 || data_type == 2
    out = 1;
elseif data_type == 3 || data_type == 4 || data_type == 5 
    out = 2;
else
    out = 3;
end

end