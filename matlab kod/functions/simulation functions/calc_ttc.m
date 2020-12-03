function ttc = calc_ttc(sa,sb,r)
% computes ttc assuming constant velocity

% extract variables form states
A0 = sa.pos;
B0 = sb.pos;
stepsize0  = [sa.speed,sb.speed];
theta0     = [sa.theta,sb.theta];

if A0<B0
    delta0 = A0-B0;
    delta_step = stepsize0(1)*exp(1i*theta0(1)) + stepsize0(2)*exp(1i*theta0(2));
    eta = real(delta0*conj(delta_step))/norm(delta_step)^2;
    %r_scaled = sqrt(norm(delta0)^2 - eta^2*norm(delta_step)^2 )/2;
    mu = (norm(delta0)^2 - 4*r^2)/norm(delta_step)^2;

    if eta^2 >= mu
        ttc = min(-eta + [-1 1]*sqrt(eta^2 - mu));%exp(r_scaled/r + 1);
    else
        ttc = inf;
    end

else
    ttc = Inf;%-eta*r_scaled*exp(r_scaled/r+1);
end

end