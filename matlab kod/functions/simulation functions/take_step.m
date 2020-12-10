function new_state = take_step(S, desired_speed, desired_angle, max_delta)
% Takes a step given current state and desired state.

% discrepency between desired and current angle:
des_delta_theta = desired_angle - S.theta;
% compute how the driver wants to steer wheel to avoid collision:
des_delta2_theta = des_delta_theta - S.dtheta;

% ensure maximum change is not violated
delta2_theta = cap(des_delta2_theta, max_delta.d2theta);
delta_theta = S.dtheta + delta2_theta;
% ensure maximum change is not violated
delta_theta = cap(delta_theta, max_delta.dtheta);

% compute new theta
theta = S.theta + delta_theta;

% desired change in speed
des_delta_speed = desired_speed - S.speed;
% ensure maximum change is not violated
delta_speed = cap(des_delta_speed, max_delta.dspeed);
% compute new speed
speed = S.speed + delta_speed;

% compute next position based on new angle and speed
A = S.pos + speed*exp(1i*theta);

% update state
[S.pos, S.theta, S.speed, S.dtheta] = deal(A,theta,speed,delta_theta);
new_state = S;

end