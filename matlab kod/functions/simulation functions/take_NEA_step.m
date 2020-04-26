function out = take_NEA_step(A0,B0, stepsize, theta, driver_prop)
%driver_prop = {stepsize0, var_step, var_theta, stability_fac, min_stepsize};
stepsize0 = driver_prop{1};
var_step = driver_prop{2};
var_theta = driver_prop{3};
stability_fac = driver_prop{4};
min_stepsize = driver_prop{5};
theta_mod_par = driver_prop{6};
speed_mod_par = driver_prop{7};

d1   = theta_mod_par(1);
p1   = theta_mod_par(2);
amp1 = theta_mod_par(3);
d2   = speed_mod_par(1);
p2   = speed_mod_par(2);
amp2 = speed_mod_par(3);

%delta_theta = amp1*s_function(real(A0),d1,p1) ^ (amp2*s_function(theta(1), d2,p2));
expected_theta = s_function(real(A0),d1,p1,amp1);
expected_speed = stepsize0(1) - amp2*exp(-p2*(theta(1)-d2)^2);

%theta(1) = theta(1) + delta_theta + normrnd(-stability_fac*(theta(1) - delta_theta), var_theta);
theta(1) = theta(1) + normrnd(-stability_fac*(theta(1) - expected_theta),var_theta);
theta(2) = theta(2) + normrnd(-stability_fac*theta(2),var_theta);
stepsize(1) = stepsize(1) + normrnd(  -stability_fac*(stepsize(1) - expected_speed),     var_step); 
stepsize(2) = stepsize(2) + normrnd(0,var_step);
stepsize = max(min_stepsize, stepsize);
A1 = A0 + stepsize(1)*exp(1i*theta(1));
B1 = B0 - stepsize(2)*exp(1i*theta(2));
out = {A1,B1,stepsize,theta};
end