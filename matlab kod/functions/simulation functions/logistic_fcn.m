function out = logistic_fcn(x,location,scale,limits)
% logistic function. 
% Goes from 0 to 1 as x tends to infty if sign is (+)
% Goes from 1 to 0 as x tends to infty if sign is (-)

% set default arguments
if nargin == 1
    scale = 1;
    limits = [0,1];
elseif nargin == 2
    location = 0;
    limits = [0,1];
elseif nargin == 3
    limits = [0,1];
end

out = (limits(2)-limits(1))./(1 + exp((x - location)/scale) ) + limits(1);

end