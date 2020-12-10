function out = logistic_fcn(x,location,scale)
% logistic function. 
% Goes from 0 to 1 as x tends to infty if sign is (+)
% Goes from 1 to 0 as x tends to infty if sign is (-)

% set default arguments
if nargin == 2
    scale = 1;
elseif nargin == 1
    location = 0;
    scale = 1;
end
    

out = 1 ./(1 + exp((x - location)/scale) );

end