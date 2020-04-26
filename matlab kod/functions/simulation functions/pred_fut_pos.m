function [Apred,Bpred] = pred_fut_pos(A0, B0,stepsize,theta)
% function that predicts the future positions of each driver using the
% constant movement extrapolation hypothesis.
Apred = A0 + stepsize(1)*exp(1i*theta(1));
Bpred = B0 - stepsize(2)*exp(1i*theta(2));
end
