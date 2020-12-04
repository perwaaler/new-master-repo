function angle = expected_angle_A(S, x0, x1)
x = real(S(1).pos);
if x<x0
    angle = 0;
elseif x>x1
    angle = pi/2;
else
    angle = asin((x-x0)/(x1-x0));
end
end