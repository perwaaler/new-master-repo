function out = calc_speed_angle(state)
% takes a 4 by 4 matrix with previous and current positions
% and returns speed and angle. row 1 corresponds to A, 
% row 2 corresponds to B.
Aold = state(1,1);
A0   = state(1,2);
Bold = state(2,1);
B0   = state(2,2);

velocityA = A0 - Aold;
velocityB = B0 - Bold;

out(1).speed = norm(velocityA);
out(2).speed = norm(velocityB);

out(1).theta = angle(velocityA);
out(2).theta = angle(velocityB);

end