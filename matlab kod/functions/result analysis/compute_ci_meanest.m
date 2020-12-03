function ci = compute_ci_meanest(data_matrix, N)
% Computes 95% confidence intervals for estimate of
% mean value m using sample data_matrix. If data_matrix is a matrix with
% each column corresponding to a sample, then function returns a matrix
% where row i is the ci correspoding to sample i. N is sizer of each
% sample.

ci = mean(data_matrix)' +  sqrt( var(data_matrix)/N )'*1.96*[-1 1];
end

