function new_state = momentum_step(S)
% computes next position based on current position, speed and
% angle. Shows where driver is in next frame if no changes to behaviour are
% made, i.e. steering wheel and amount of braking are held constant.

RUprop = S.RUprop;

% speed prediction
S.speed  = max(S.speed + S.dspeed, 0);
% add prediction error and ensure non-negative
if S.speed > 0
    S.speed  = gam_rnd(S.speed, RUprop.pred_speed_err(S.id));
end

% theta prediction
S.theta  = S.theta + S.dtheta;
% add prediction error
S.theta = normrnd(S.theta, RUprop.pred_theta_err(S.id));
S.pos    = S.pos   + S.speed*exp( 1i*S.theta );

new_state = S;

end