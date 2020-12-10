function ttc = calc_ttc(S,r)
% computes ttc assuming constant velocity

% extract variables from state
if isa(S,'struct') %input is a structure with fields
        A0 = S(1).pos;
        B0 = S(2).pos;
        speed0 = [S(1).speed,S(2).speed];
        theta0 = [S(1).theta,S(2).theta];
        
elseif isa(S,'double') %input is a 4by4 matrix with old and current positions
    A0 = S(1,2);
    B0 = S(2,2);
    velocity = calc_speed_angle(S);
    speed0   = [velocity(1).speed,velocity(2).speed];
    theta0   = [velocity(1).theta,velocity(2).theta];
end

% calculate TTC
deltaP = A0-B0;
deltaV = speed0(1)*exp(1i*theta0(1)) - speed0(2)*exp(1i*theta0(2));
eta = real(deltaP*conj(deltaV))/norm(deltaV)^2;
mu = (norm(deltaP)^2 - 4*r^2)/norm(deltaV)^2;

if eta^2 >= mu
    ttc = min(-eta + [-1 1]*sqrt(eta^2 - mu));
else
    ttc = inf;
end



end