function phi = orient2action(orient)
% returs a character string saying how the rad users should steer in order
% to increase t_min.

phi = ["1","2"];
for i=1:2
    if orient(i)=="left"
        phi(i) = "turn right";
    else
        phi(i) = "turn left";
    end
end
end
