function cost = cost_fcn(x,d,p)
% function used to calculate the cost associated with selection of 
% a particular threshold. See thr_autofind.
cost = zeros(1,length(x));
cost(abs(x)<d/2) = abs(   x( abs(x)<d/2 ) /d);
cost(abs(x)>=d/2) = 1/2 - 1/2^p + abs( x( abs(x)>=d/2 )/d ).^p;
end