function plotbeta(interval,mean,sigma,res)
% plot beta pdf
if nargin == 3
    res = 400;
end
x = linspace(interval(1),interval(2),res);
[a,b] = betapar(mean,sigma);
plot(x, betapdf(x, a, b))
end