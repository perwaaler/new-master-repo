function [A0,B0] = calc_prev_pos(S)
% computes the previous poitions based on current state
[A,B] = S.pos;
[speedA,speedB] = S.speed;
[thetaA,thetaB] = S.theta;

A0 = A - speedA*exp(1i*thetaA);
B0 = B - speedB*exp(1i*thetaB);
end
