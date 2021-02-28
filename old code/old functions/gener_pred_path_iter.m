function pos = gener_pred_path_iter(state,max_delta,k)
% iteratively computes the k next positions assuming constant 
% change of speed and angle.

pos = zeros(1,k+1);
pos(1) = state.pos;


if state.dspeed<=0
    % k0 is the maximum number of steps RU can take before speed becomes
    % negative
    ks = floor(abs(state.speed/state.dspeed));
    ks = min(k,ks);
else
    ks = k;
end

speed = state.speed + (1:ks)*state.dspeed;

kt = floor(abs(pi*4/5/state.dtheta));
theta = state.theta + (1:kt)*state.dtheta;
theta = [theta, ...
        (state.theta + sign(state.dtheta)*pi*3/4)*ones(1,k-kt)];

for i=1:ks
    pos(i+1) = pos(i) + speed(i)*...
               exp(1i*theta(i));
end
pos = pos(1:ks);

end