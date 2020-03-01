function [param, p_u] = bootstrap_est(data, Nenc, init_val, u, shake_var)
% function that generates bootstrapped sample and estimates the parameters
% for the bootstrapped sample

resampling = randsample(Nenc, Nenc, true);              % indeces used to bootstrap
data_bs = data(resampling,:);                           % bootstrapped sample
p_u = sum(data_bs(:)>u)/length(data_bs(:));

param = est_par_w_shake( data_bs, init_val, u, shake_var);
if param == init_val
    param = [nan nan];
end

end