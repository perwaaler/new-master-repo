function out = lognormal_par(mu,sigma)
% function that takes the desired mean and standard deviation of the 
% lognormal distribution and transforms them into parameters to plug into
% the lognormal function, ex: lognrnd(mu,sigma).
mu_ln = log(mu^2/sqrt(mu^2+sigma^2));
sigma_ln = sqrt( log(1 + sigma^2/mu^2) );
out = [mu_ln, sigma_ln];
end