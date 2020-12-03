function param = est_par(data, init_val, u)
% function that estimates the parameter of the GPD distribution of the
% exceedences over a threshold u

% find exceedences over threshold u
exceed = data(data > u); 
% negative loglikelihood for the data
negL = @(par) -sum( log( gppdf(exceed, par(2), par(1), u)) );
%find the optimal parameters
param = fminsearch(negL, init_val);

end