function line_coeff = path_tangent(state)
% function that takes the state S and computes the coefficients of
% the tangent line of the trajectory at the position of S. Also returns the
% x-component of the gradient vector.
pos_prev = state.pos - state.speed*exp(1i*state.theta);
pos_curr = state.pos;

diff = pos_curr - pos_prev;
x_gradient = real(diff);
slope = imag(diff)/real(diff);
intercept = imag(pos_curr) - slope*real(pos_curr);
% y1 = slope*x1 + c --> c = y1 - k*x1
line_coeff = [intercept,slope,x_gradient];

end