function [p,t_min,d_min] = p_engage_EA_v2(A0,B0,stepsize,theta, driver_prop, detec_par, farsight)
% function that computes the probability of EA mode being engaged
stepsize0 = driver_prop{1};
stability_fac = driver_prop{4};
theta_mod_par = driver_prop{6};
speed_mod_par = driver_prop{7};

d1   = theta_mod_par(1);
p1   = theta_mod_par(2);
amp1 = theta_mod_par(3);
d2   = speed_mod_par(1);
p2   = speed_mod_par(2);
amp2 = speed_mod_par(3);

shift_d = detec_par(1);
pow_d   = detec_par(2);
shift_t = detec_par(3);
pow_t   = detec_par(4);
amp     = detec_par(5);
%future_pos = zeros(6,2);
%future_pos(1,:) = [A0,B0];
dist = zeros(1,farsight+1);
dist(1) = norm(A0-B0);

for k=1:farsight
    expected_theta = s_function(real(A0),d1,p1,amp1);
    expected_speed = stepsize0(1) - amp2*exp(-p2*(theta(1)-d2)^2)*(s_function(1,-0.6,2)+1);
    
    theta(1) = theta(1) - stability_fac*( theta(1) - expected_theta );
    theta(2) = theta(2) - stability_fac*theta(2);
    
    stepsize(1) = stepsize(1) - stability_fac*(stepsize(1) - expected_speed);
   
    [A0,B0] = pred_fut_pos(A0, B0, stepsize, theta);
    dist(1+k) = norm(A0-B0);
 %   plot_pos(A0, B0, 0.01, 6, 0.3, 1) 
 %   future_pos(k,:) = [A0,B0];
end
[d_min, t_min] = min(dist);

p = 1./(1+exp(-(2 - 0.17*t_min - 1.9624*d_min)));
end
