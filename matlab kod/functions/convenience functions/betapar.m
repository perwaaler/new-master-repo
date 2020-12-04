function [alfa,beta] = betapar(meanx,sigmax)
mu = (meanx.*(1 - meanx)./sigmax.^2 - 1);
alfa = meanx.*mu;
beta = (1 - meanx).*mu;
end