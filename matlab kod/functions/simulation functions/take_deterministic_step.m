function S = take_deterministic_step(S)
% function that computes the step nex step using no randomization.
S.theta = cap(S.theta + S.dtheta);
S.speed = S.speed + S.dspeed;
S.pos = S.pos + S.speed*exp(1i*S.theta);
end