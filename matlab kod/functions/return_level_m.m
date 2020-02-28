function x_m = return_level_m(m, sigma, xi, u, p_exceed)
% function that computes the mth obs. return level, which is the level x_m
% that on average gets exceeded once every m observations.
x_m = u + sigma./xi.*( (m*p_exceed).^xi - 1); 
end