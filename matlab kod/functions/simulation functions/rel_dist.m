function out = rel_dist(S)
% calculates relative distance
A_old = S(1).pos - S(1).speed*exp(1i*S(1).theta);
B_old = S(2).pos - S(2).speed*exp(1i*S(2).theta);

d_old = norm(A_old - B_old);
d = norm(S(1).pos - S(2).pos);
out = d/d_old;
end
