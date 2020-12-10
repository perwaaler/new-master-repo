function out = gam_rnd(meanx,sigmax,n_rnd)
% same as gamrnd, but takes mean and standard deviation as parameters

if nargin == 2
    n_rnd = 1;
end

% reparametrize
k = meanx.^2./sigmax.^2;
theta = sigmax.^2./meanx;

out = gamrnd(k,theta,1,n_rnd);

end