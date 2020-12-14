function plotgamma(interval,mean,sigma,res)

if nargin == 3
    res = 300;
end

x = linspace(interval(1),interval(2),res);
[a,b] = gammapar(mean,sigma);
plot(x, gampdf(x, a, b))
end