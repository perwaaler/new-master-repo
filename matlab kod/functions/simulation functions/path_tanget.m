function line_coeff = path_tanget(state)
% function that takes the state S and computes the coefficients of
% the tangent line of the trajectory at the position of S.
pos_prev = state.pos - state.speed*exp(1i*state.theta);
pos_curr = state.pos;

diff = pos_curr - pos_prev;
slope = imag(diff)/real(diff);
intercept = imag(pos_curr) - slope*real(pos_curr);
% y1 = slope*x1 + c --> c = y1 - k*x1
line_coeff = [intercept,slope];

end