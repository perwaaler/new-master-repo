function tstop = predict_tstop(state)
% takes the current state and predicts the amount of time untill Road user
% will stop moving, rounded down to the nearest integer.
if state.dspeed>0
    tstop = inf;
elseif state.dspeed<=0
    if state.speed==0
        tstop = 0;
    else
        tstop = floor(abs(state.speed/state.dspeed));
    end
end
end