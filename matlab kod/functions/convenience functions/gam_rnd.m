function out = gam_rnd(meanx,sigmax,n_rnd)
% same as gamrnd, but takes mean and standard deviation as parameters

if nargin == 2
    n_rnd = 1;
end

[k,theta] = gammapar(meanx,sigmax);
out = gamrnd(k,theta,1,n_rnd);

end