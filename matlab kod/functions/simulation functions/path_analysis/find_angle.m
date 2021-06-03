function angle_t = find_angle(paths,Id,t,pp,h)
% function that returns the angle of the path at time t, using
% extrapolation based on the piecewise polynomial fit pp.
% * paths is a field containing information about generated paths
% * Id is the road users id
% * t is the time at which the slope is to be computed
% * pp is the piecewise polynomial fitted to the path
% * h is the x-increment used to approimate the derivative

der(:,1) = differentiate(pp, t,-h);
der(:,2) = differentiate(pp, t,h);


[m,I] = max(sum(abs(der)));

if m==0
    % this condition implies that the road user is standing still at time t
    angle_t = paths.theta(Id);
else
    % road user is moving at time t
    angle_t = angle( der(1,I) + 1i*der(2,I) );
end

end