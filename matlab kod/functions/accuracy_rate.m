function acc_rate = accuracy_rate(estimate, true_val, cut_off)
% Estimates the probability that the estimator deviates from the true value
% by less that cut_off*100 percent.

est_to_truth_ratio = estimate/true_val;
acc_rate = mean( abs(est_to_truth_ratio - 1) < cut_off);

end