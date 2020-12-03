function out = logpar(mux,sigmax)
% function that takes desired mean and standard deviation and returns the
% corresponding lognormal parameters
mu =    log( mux^2/sqrt(mux^2 + sigmax^2) );
sigma = log( 1 + sigmax^2/mux^2 );
out = [mu, sigma];
end