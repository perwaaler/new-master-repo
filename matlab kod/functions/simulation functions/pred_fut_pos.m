function [Apred,Bpred] = pred_fut_pos(A0, B0,stepsize,theta0)
Apred = A0 + stepsize(1)*exp(-1i*theta0(1));
Bpred = B0 - stepsize(2)*exp(-1i*theta0(2));
end
