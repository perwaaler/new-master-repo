function plot_logistic_fcn(interval,location,scale,res)
% plots the logistic function over desired interval

% set default parameters
if nargin==1
    location = 0;
    scale = 1;
    res = 300;
    
elseif nargin==3
    res = 300;
end

x = linspace(interval(1),interval(2),res);
plot(x, logistic_fcn(x, location, scale))
end