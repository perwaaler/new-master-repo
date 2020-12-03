function new_state = take_step(A0, speed0, angle0, dtheta0, desired_speed, desired_angle, max_change)
max_delta_speed = max_change(1);
max_delta_theta = max_change(2);
max_delta_dtheta = deg2rad(4);

angle_diff = desired_angle - angle0;
%angle_acc  = angle_diff - dtheta0;
speed_diff = desired_speed - speed0;

angle_diff_diff =  sign(angle_diff)*min(max_delta_dtheta, abs(angle_diff - dtheta0));


angle_diff = dtheta0 + angle_diff_diff;

% compute new angle and speed, and take next step
angle1 = angle0 + sign(angle_diff)*min(max_delta_theta, abs(angle_diff));
speed1 = speed0 + sign(speed_diff)*min(max_delta_speed, abs(speed_diff));
A1 = A0 + speed1*exp(1i*angle1);
new_state = [A1, speed1, angle1];
end