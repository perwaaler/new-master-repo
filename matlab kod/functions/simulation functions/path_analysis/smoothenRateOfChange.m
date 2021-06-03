function [S1,S2] = smoothenRateOfChange(enc,stat,ReachBack)

% how far back in time to take average
count_back = max(1,stat.frame - ReachBack);
% extract states
prev_states = enc.states(count_back:stat.frame); 
% 
nSmoothInd = size(prev_states,1);

dtheta = [0,0];
dspeed = [0,0];
for i=1:nSmoothInd
    dtheta(1) = dtheta(1) + prev_states{i}(1).dtheta;
    dtheta(2) = dtheta(2) + prev_states{i}(2).dtheta;
    dspeed(1) = dspeed(1) + prev_states{i}(1).dspeed;
    dspeed(2) = dspeed(2) + prev_states{i}(2).dspeed;
end
dtheta = dtheta/nSmoothInd;
dspeed = dspeed/nSmoothInd;

S1 = prev_states{end}(1);
S2 = prev_states{end}(2);

S1.dspeed = dspeed(1);
S2.dspeed = dspeed(2);

S1.dtheta = dtheta(1);
S2.dtheta = dtheta(2);

end