function [ci_sigma_u, ci_xi_u, ci_pc_u] = compute_bs_ci(param, pc, sigma_sample, xi_sample, pc_sample)

sigma_sample(isnan(sigma_sample)) = [];
xi_sample(isnan(xi_sample)) = [];
pc_sample(isnan(pc_sample)) = [];

se_sigma = std(sigma_sample);
se_xi = std(xi_sample);
se_p_nea = std(pc_sample);

ci_sigma = sort(  param(1) + [-1 1]*1.96*se_sigma  );
ci_xi =    sort(  param(2) + [-1 1]*1.96*se_xi     );
ci_p_nea = sort(  pc + [-1 1]*1.96*se_p_nea     );

ci_sigma_u = ci_sigma';
ci_xi_u = ci_xi';
ci_pc_u = ci_p_nea';

end