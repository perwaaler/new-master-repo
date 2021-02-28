
alpha = 0.5;
m=100;
kkk = 1000;
x = linspace(0,50,kkk);
sigma1 = 1.8;
xi1 = .08;
sigma2 = 1.5;
xi2 = .15;

clf
plot(x, gpcdf(x,xi1, sigma1),'b')
hold on
plot(x, gpcdf(x,xi2,sigma2),'r')
plot(x, gpcdf(x,xi2,sigma2) - gpcdf(x,xi1,sigma1),'r')
plot(x,zeros(1,kkk))

return_level_m(m, sigma1, xi1, 0, 1)
return_level_m(m, sigma2, xi2, 0, 1)

p01 = 1-gpcdf(30,xi1,sigma1)
p02 = 1-gpcdf(30,xi2,sigma2)
%%
clf
for i=1:4
    hold on
    subplot(211)
    hold on
    plot(x, plot(x, gpcdf(x,xi2^i,sigma2^2),'r'))
    subplot(212)
    hold on
    plot(x, plot(x, gpcdf(x,xi2^i,sigma2),'r'))
    pause(1)
end
%%
s=tf('s');
z=[3 6 12];
subplot(2,1,1)
hold on
for i=1:3
    zc=z(i);
    g=(15/zc)*((s+zc)/(s^2+3*s+15));
    step(g)
end
legend('z = 3','z = 6','z = 12')
subplot(2,1,2)
hold on
for i=1:3
    zc=z(i);
    g=(15/zc)*((s+zc)/(s^2+3*s+15));
    impulse(g)
end
legend('z = 3','z = 6','z = 12')
