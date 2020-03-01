function [A1, B1, theta] = take_evasive_step(A0,B0,step_par,stepsize,var_step, min_stepsize,theta0,thetavar,stability_fac,thetamod_k,thetamod_p,thetamod_s)
% take next step
pred_posA = A0 + step_par(1)*step_par(2); % prediction of the A's future position
pred_posB = B0 - step_par(1)*step_par(2); % prediction of the B's future position
pred_diff = pred_posA - pred_posB; % difference vector between predicted future position of A and B
pred_dist = norm(pred_diff);       % distance between predicted positions at time i+1
stepsize = stepsize + normrnd([0,0],var_step);
stepsize(1) = max(min_stepsize,stepsize(1));
stepsize(2) = max(min_stepsize,stepsize(2));
theta = normrnd(theta0,thetavar) + normrnd(-stability_fac*theta0,thetavar) + exp(-0.2*pred_dist)*thetamod_k*sign(imag(pred_diff))*exppdf(abs(imag(thetamod_s*pred_diff))^thetamod_p, 1);
A1 = A0 + stepsize(1)*exp(1i*theta(1));
B1 = B0 - stepsize(2)*exp(1i*theta(2));
end