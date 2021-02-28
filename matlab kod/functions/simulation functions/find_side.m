function out = find_side(u,v)
% function that figures out on which side of the vector u that vector v is
% on, left or right. u and v are given as complex numbers.
% u can also be given by the a vector

phi   = mod(angle(u),2*pi); % rad2deg(phi)
theta = mod(angle(v),2*pi); % rad2deg(theta) rad2deg(phi-pi)

if phi<pi
    
    if phi < theta &&  theta < phi + pi
        out = "left";
    else
        out = "right";
    end
    
elseif phi>pi
    
    if phi-pi > theta ||  phi < theta
        out = "left";
    else
        out = "right";
    end
end

end