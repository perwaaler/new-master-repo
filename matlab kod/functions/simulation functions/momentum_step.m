function new_state = momentum_step(S,max_delta,free_theta)
% computes next position based on current position, speed and
% angle. Shows where driver is in next frame if no changes to behaviour are
% made, i.e. steering wheel and amount of braking are held constant.

RUprop = S.RUprop;
    
if nargin==2
    % allow theta to vary
    free_theta=1;
end

% ensure non-negative speed
S.speed  = max(S.speed + S.dspeed, 0);

% theta prediction
S.dtheta = cap(S.dtheta, max_delta.dtheta*S.speed^.8);
S.theta  = S.theta + S.dtheta*free_theta;
% add prediction error
S.theta = normrnd(S.theta, RUprop.pred_theta_err(S.id));
S.pos    = S.pos   + S.speed*exp( 1i*S.theta );

new_state = S;


end