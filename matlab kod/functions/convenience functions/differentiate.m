function out = differentiate(f,x,h)
% differentiates the function f at the point x
out = (f(x+h) - f(x))/h;
end