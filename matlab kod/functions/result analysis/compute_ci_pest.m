function ci = compute_ci_pest(prob_est, N)
% computes 95% confidence intervals for estimate of a
% probability p using sample of size N. If prob_set is a vector, then it
% returns a k by 2 matrix with row i corresponding to estimate i.

if size(prob_est,1)<size(prob_est,2)
    prob_est=prob_est';
end
ci = prob_est +  sqrt( prob_est.*(1-prob_est)/N )*1.96*[-1 1];

end

