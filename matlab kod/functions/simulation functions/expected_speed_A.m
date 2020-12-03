function speed = expected_speed_A(x, x0, x1, init_speed)
if x<x0
    speed = init_speed;
elseif x>x1
    speed = init_speed;
else
    %speed = init_speed + 0.2*init_speed*(cos(4*asin((x-x0)/(x1-x0)))-1);
%     angle = 4*asin((x-x0)/(x1-x0));
    speed = init_speed + 0.02*init_speed*(x-x0)*(x-x1);
end
end