function out = differentiate(f,x,h)
% Approximates derivative of the function f at the point x.
out = (f(x+h) - f(x))/h;
end