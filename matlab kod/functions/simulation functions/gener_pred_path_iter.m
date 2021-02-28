function path = gener_pred_path_iter(state,max_delta,k,sparse_ind)
% iteratively computes the k next positions assuming constant 
% change of speed and angle. tstop is the predicted time untill RU is
% standing still. If he is accelarating this time is inf.

pos = zeros(1,k+1);
pos(1) = state.pos;
theta0 = state.theta;
% variable that prevents the path from looping
free_theta = 1;

if state.dspeed>0
    tstop = inf;
elseif state.dspeed<=0
    if state.speed==0
        tstop = 0;
    else
        tstop = floor(abs(state.speed/state.dspeed));
    end
end

for i=1:k
    if abs(state.theta-theta0)>=2.6180
        % force theta to be fixed after turning more than 2.6
        free_theta = 0;
    end
    
    state = momentum_step(state,max_delta(1),free_theta);
    
    if state.speed==0
        % RU has stopped. return length of path
        k = i;
        break
    else
        pos(i+1) = state.pos;
    end
    
end

sparse_ind = sparse_ind(sparse_ind<=k);

path.pos = pos(sparse_ind);
path.ind = sparse_ind;
path.slope = tan(state.theta);
path.length = k;
path.tstop = tstop;

% startslope = complex2cart(exp(1i*theta0));
% endslope   = complex2cart(exp(1i*state.theta));
% slopes = [startslope, endslope,times];

end