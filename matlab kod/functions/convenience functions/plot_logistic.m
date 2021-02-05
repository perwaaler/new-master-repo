function plot_logistic(interval,location,scale,limits,res)
% plots the logistic function over desired interval

% set default parameters
if nargin==1
    location = 0;
    scale = 1;
    limits = [0,1];
	res = 300;
elseif nargin==2
    scale = 1;
    limits = [0,1];
    res = 300;
elseif nargin==3
    limits = [0,1];
    res = 300;
elseif nargin==4
    res = 300;
end

x = linspace(interval(1),interval(2),res);
plot(x, logistic_fcn(x, location, scale, limits))
end