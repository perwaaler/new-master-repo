function side_lines = find_path_boundaries(state,r)
% function that finds the parameters for the lines that defines the path of
% the local trajectory.

% find parameters of line:
A0 = state.pos - state.speed*exp(1i*state.theta);
A1 = state.pos;

x0 = real(A0);
y0 = imag(A0);

% find parameters of the tangent line of the trajectory at current frame
line_par= find_line_par(A0,A1);
alfa = line_par(2);

deltax = r*alfa/sqrt(alfa^2 + 1);
deltay = deltax/alfa;

% find points that line1 and line2 must pass through
ux = x0 + deltax;
uy = y0 - deltay;
vx = x0 - deltax;
vy = y0 + deltay;

beta_u = uy - alfa*ux;
beta_v = vy - alfa*vx;

line_u = [beta_u,alfa];
line_v = [beta_v,alfa];

side_lines = [line_u;line_v];

% plot_circle(A0,r)
% hold on
% drawline(line_par)
% drawline([beta_u,alfa])
% drawline([beta_v,alfa])
% hold on
% plot_circle(ux+1i*uy,r)
% plot_circle(vx+1i*vy,r)

end
