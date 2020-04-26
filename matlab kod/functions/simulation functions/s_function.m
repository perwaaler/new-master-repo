function out = s_function(x, location, power, amp)
% s function, used to construct p_EA. 
% Goes from 0 to 1 as x tends to infty if sign is (+)
% Goes from 1 to 0 as x tends to infty if sign is (-)
if nargin == 3
    out = 1  /(1 + exp(power*(x - location)) );
else
    out = amp/(1 + exp(power*(x - location)) );
end
end