function angle = E_theta_A(S, RUprop)
% computes the expected direction of A as a function of x

x  = real(S(1).pos);
x0 = RUprop.turnstart;
x1 = RUprop.turnend;

if x<x0
    angle = 0;
elseif x>x1
    angle = pi/2;
else
    angle = asin((x-x0)/(x1-x0));
end

end