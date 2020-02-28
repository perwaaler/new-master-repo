function param = est_par(data, init_val, u)
% function that estimates the parameter. Takes 

exceed = data(data > u);
negL = @(par) -sum( log( gppdf(exceed, par(2), par(1), u)) );
param = fminsearch(negL, init_val);

end