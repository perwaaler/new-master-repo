function speed = E_speed_A(S, RUprop)
% computes the expected speed of A as a function of x.

x = real(S(1).pos);
x = x + RUprop.start_brake;
x0        = RUprop.turnstart;
x1        = RUprop.turnend;
avg_speed = RUprop.avg_speed(1);
brake     = RUprop.turn_brake;

if x<x0
    speed = avg_speed;
elseif x>x1
    speed = avg_speed;
else
    speed = avg_speed + brake*avg_speed*(x-x0)*(x-x1);
end
end