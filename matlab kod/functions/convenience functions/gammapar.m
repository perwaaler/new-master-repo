function [k,theta] = gammapar(meanx,sigmax)
k = meanx^2/sigmax^2;
theta = sigmax^2/meanx;
end