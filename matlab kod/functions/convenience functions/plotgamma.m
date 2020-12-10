function plotgamma(interval,mean,sigma,res)
x = linspace(interval(1),interval(2),res);
[a,b] = gammapar(mean,sigma);
plot(x, gampdf(x, a, b))
end