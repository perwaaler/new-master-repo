function angle = expected_angle_A(x, x0, x1)
if x<x0
    angle = 0;
elseif x>x1
    angle = pi/2;
else
    angle = asin((x-x0)/(x1-x0));
end
end