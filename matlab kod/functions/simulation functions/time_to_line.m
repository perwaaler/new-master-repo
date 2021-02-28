function ttl = time_to_line(state,line_par,r)
% computes the time that elapses until a ball reaches the line with 
% given by line_par, assuming constant velocity.

% alfa = 1.2;
% beta = 1;
% clf
% drawline([beta,alfa])
% 
% A0 = -1+4i;
% A1 = -.5 + 2.5i;
% plot_circle(A0,r)
% plot_circle(A1,r)

beta = line_par(1);
alfa = line_par(2);

% extract and compute previous and current positions
A0 = state.pos - state.speed*exp(1i*state.theta);
A1 = state.pos;

u = real(A0);
v = imag(A0);
u1 = real(A1);
v1 = imag(A1);

du = u1-u;
dv = v1-v;

A = v - beta - alfa*u;
B = dv - alfa*du;
C = u + alfa^-1*(beta - v);
D = du - alfa^-1*dv;

phi   = A*B + C*D;
theta = B^2 + D^2;
eta   = A^2 + C^2;

% compute time to reach line
t = -phi/theta + [-1 1]*sqrt((phi/theta)^2 - ...
                             (eta - r^2*(alfa^-1+alfa)^2 )/theta);

ttl = min(t);  

% if ttl<0
%     msg = 'ball is moving away from line';
%     error(msg)
% end
% (alfa^-1+alfa)^-2*((A + t(1)*B)^2 + (C + D*t(1))^2 )
% P_t = A0 + t(1)*(A1-A0)
% plot_circle(P_t,r,'yellow')
% norm(P_t - (x+y*1i))
end

