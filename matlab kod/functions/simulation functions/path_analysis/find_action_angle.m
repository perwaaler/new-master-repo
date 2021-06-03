function phi = find_action_angle(orient)
% takes orientation and determines how to change direction in order to
% avoid collision. Phi is a vector that contais the action angles.

phi = nan(1,2);
for i=1:2
    if orient(i)=="left"
        phi(i) = -1;
    else
        phi(i) = +1;
    end
end

end