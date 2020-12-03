function describe_variable(variable)
% function that outputs a string describing the variable
if variable == 'driver_prop'
    sprintf(' stepsize0 \n var_step \n var_theta\n stability_fac\n min_stepsize\n theta_mod_par\n speed_mod_par')
elseif variable == 'state_A'
    sprintf(' A0 \n speed_A \n theta_A \n Dtheta_A')
end

