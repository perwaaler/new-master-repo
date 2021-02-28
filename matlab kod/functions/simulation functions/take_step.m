function new_state = take_step(S, desired_speed, desired_angle, max_delta)
% Takes a step given current state and desired state.

% if only S is given, predict next step using constant dtheta and dspeed
if nargin == 1
    new_state = momentum_step(S,max_delta);
    
else
% discrepency between desired and current angle:
des_delta_theta = desired_angle - S.theta;
% compute how the driver wants to steer wheel to approach desired angle:
des_delta2_theta = des_delta_theta - S.dtheta;

% ensure maximum change is not violated
delta2_theta = normrnd(cap(des_delta2_theta, max_delta.d2theta),...
                                            S.RUprop.d2theta_std(S.id));
delta_theta = S.dtheta + delta2_theta;
% ensure maximum change is not violated
delta_theta = cap(delta_theta, max_delta.dtheta*S.speed^.8);

% compute new theta
theta = S.theta + delta_theta;

% desired change in speed
des_delta_speed = desired_speed - S.speed;
% ensure maximum change is not violated
delta_speed = cap(des_delta_speed, max_delta.dspeed);
% compute new speed
speed = max(0,S.speed + delta_speed);
delta_speed = S.speed - speed;

% compute next position based on new angle and speed
A = S.pos + speed*exp(1i*theta);

% update state
S.pos    = A;
S.theta  = theta; 
S.speed  = speed;
S.dtheta = delta_theta;
S.dspeed = delta_speed;
new_state = S;

end

end