function out = pw_fcn(pp,x,tend)
% function used in instances when a road user stops. Used specifically for
% optimization in fminsearch to find t_min or proximity times. 
% tend is the amount of time untill RU is predicted to stop.
% pp is a piecewise polynomial fitted to the predicted future positions.
if x<=0
    out = pp(0);
elseif 0<x && x<=tend
    out = pp(x);
else
    out = pp(tend);
end
end
    
