function pred_path = gener_pred_path(state,k)
% generates a path that assumes constant change of direction and speed.
% Takes k steps

% 
% A0 = 1 + 2i;
% theta0 = 0.1;
% dtheta = deg2rad(5);
% speed0 = 0.16;
% dspeed = -0.003;
% plot_circle(A0,0.3)
% k = 100;


A0 = state.pos;
theta0 = state.theta;
dtheta = state.dtheta;
speed0 = state.speed;
dspeed = state.dspeed;

% ensure non-negative speed
if dspeed<=0
    % k0 is the maximum number of steps RU can take before speed becomes
    % negative
    k0 = floor(abs(speed0/dspeed));
    k0 = min(k,k0);
    kdiff = k-k0;
else
    k0 = k;
    kdiff = 0;
end

k_step  = 1:k0;
k_one = ones(1,k0);

% lower diagonal matrix of ones
tri_one = repelem(k_one,k0,1) ;
tri_one = tril(tri_one);
tri_one = sparse(tri_one);

% lower diagonal matrix of increasing values
step_triang = repelem(k_step,k0,1);
step_triang = tril(step_triang);
step_triang = sparse(step_triang);

% create direction matrix
T1 = repelem(exp(1i*k_step*dtheta),k0,1) ;
T1 = tril(T1);
T1 = sparse(T1);

T2 = repelem(exp(1i*theta0)*k_one,k0,1) ;
T2 = tril(T2);
T2 = sparse(T2);

% create matrices containging speeds and directions
SPEED = step_triang*dspeed + tri_one*speed0;
THETA = T1.*T2;

% matrix containg information about direction at each step
DELTA = sum(SPEED.*THETA,2);

% predicted future positions:
pred_path = A0 + DELTA;

pred_path(k0:k0+kdiff) = pred_path(end)*ones(1,kdiff+1);

% clf
% plot_circle(A0,0.3)
% for i=1:k
% plot_circle(pred_path(i),0.3)
% end

end
