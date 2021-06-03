function paths = gener_pred_path_iter_test(stateA,stateB,k,sparse_ind,...
                                            genInfo)
% iteratively computes the k next positions assuming constant 
% change of speed and angle. 
% * tstop is the predicted time untill RU is standing still
% If he is accelarating this time is inf
% * stateA and stateB are initial states
% * Relevance index correpsonds to the part where there is movement, or 
% the path generated before termination criteria is met, whichever comes 
% first.

if nargin == 4
    genInfo.genPath = false;
    genInfo.plotRes = false;
    genInfo.analyze = false;
end
analyze = genInfo.analyze;

if genInfo.genPath
    states = gen_random_states();
    stateA = states{1};
    stateB = states{2};
end

%     states = gen_random_states();
%     stateA = states{1};
%     stateB = states{2};
r = stateA.RUprop.r;
posA = zeros(1,k+1);
posB = zeros(1,k+1);
posA(1) = stateA.pos;
posB(1) = stateB.pos;

thetaA = stateA.theta;
thetaB = stateB.theta;

% indicator and index variables used to determine if loop should terminate
contact         = false;
tmin_found      = false;
tprox_found     = false;
end_criter_met  = false;
tick.div        = 0; % check for divergene
tick.t_prox     = 0; % check if t_prox has been found
tick.t_min      = 0; % check if tmin has been found
% have they taken first zero step?
first_zero_step  = [0,0];
% relevance_ind saves the index every time a non-zero step is taken
relevance_ind    = [k,k];
% initiate candidates information
candidates.index = {1,1};
candidates.prime_cand = [1,1];
free_theta = [1,1];
i_min = 1;
i_first_cont = [];

dist.min       = [inf, abs(posA(1) - posB(1))];
dist.prox(1)   = dist.min(1);
dist.thr       = 2*r;
dmin           = dist.min(2);

% % vector with the times they are predicted to stop
% tstop = [predict_tstop(stateA), predict_tstop(stateB)];

for i=1:k
    
    % prevent going in circles
    free_theta = straighten_out(stateA,stateB,thetaA,thetaB,free_theta);
    
    if not(tprox_found)
        [candidates,done_intersecting] = est_prox_points(candidates,...
                                                  posA, posB,i,4,dist.thr);
        dist.prox(2) = candidates.min_dist(1);            
    end
    
    if analyze
            plot_circle(posA(i),0.3,'green');
            hold on
            plot_circle(posB(i),0.3,'blue');
        if tprox_found
            plot_circle(posA(candidates.prime_cand(1)),.3,'black',2);
            plot_circle(posB(candidates.prime_cand(2)),.3,'black',2);
        end
        if tmin_found
            plot_circle(posA(i_min),.3,'magenta',2);
            plot_circle(posB(i_min),.3,'magenta',2);
        end
    end

    %%%% UPDATE STATES AND TAKE NEW STEP %%%%
    if first_zero_step(1)==0
        % first zero step has not been taken, take new step
        stateA = momentum_step(stateA, free_theta(1));
        posA(i+1) = stateA.pos;
    else
        % road user is standing still, do not update
        posA(i+1) = stateA.pos;
    end
    if first_zero_step(2)==0
        % first zero step has not been taken, take new step
        stateB = momentum_step(stateB, free_theta(2));
        posB(i+1) = stateB.pos;
    else
        % road user is standing still, do not take new step
        posB(i+1) = stateB.pos;
    end
           
    %%%% HAS THE FIRST ZERO STEP BEEN TAKEN? %%%%
    if stateA.speed==0 && first_zero_step(1)==0
        % save index that marks the end of movement
        relevance_ind(1)   = i;
        % mark that first zero step has been taken
        first_zero_step(1) = 1;
    end
    if stateB.speed==0 && first_zero_step(2)==0
        % save index that marks the end of movement
        relevance_ind(2)   = i;
        % mark that first zero step has been taken
        first_zero_step(2) = 1;
    end
    
    % update minimum distances
    if not(tmin_found)
        dist.min(1) = dist.min(2);
        dist.min(2) = abs(posB(i+1) - posA(i+1));
        
        if dist.min(2)<dmin
            i_min = i+1;
            dmin = dist.min(2);
            if dmin<dist.thr
                tmin_found = true;
            end
        end
    end
    
    %%%% LOOK FOR MOMENT OF FIRST CONTACT %%%%
    if dist.prox(2) < dist.thr && not(contact)
        % mark that the paths have overlapped
        contact = true;
        % save index for first contact
        i_first_cont = candidates.prime_cand;
    end

    %%%% ARE THE ROAD USERS GETTING CLOSER? %%%%
    if not(tmin_found) && dist.min(1) < dist.min(2)
        tick.t_min = tick.t_min + 1;
        % save the approximation of t_min
        if tick.t_min==4 || dist.min(2)<dist.thr
            tmin_found = true;
        end
    end
    
    % check if end-criteria have been met
    [end_criter_met,tick] = stop_gen_path(dist, contact,...
                        first_zero_step, done_intersecting,...
                        tick, end_criter_met, analyze);
    
    if not(tprox_found)
        [tprox_found,tick,dist] = check_if_tprox_found(tick,dist);
        dist.prox(1) = dist.prox(2);
    end
   
    
    if end_criter_met
        relevance_ind = min(i+1,relevance_ind); 
        if analyze
            plot_circle(posA(candidates.prime_cand(1)),.3,'black',4);
            plot_circle(posB(candidates.prime_cand(2)),.3,'black',4);
            plot_circle(posA(i_min),.3,'magenta',2);
            plot_circle(posB(i_min),.3,'magenta',2);
        end
        break
    end

end

% t_prox = candidates.prime_cand;
% est_prox_points(candidates,...
%                           posA, posB,i,4,dist.thr);
% pA = posA(t_prox - 1)
% pB = posB(t_prox(2))
% est_prox_points(candidates,...
%                           posA, posB,i,4,dist.thr);



% full path, including zero-speed positions. Used primarily for plotting
pathA.full = posA(1:i+1);
pathB.full = posB(1:i+1);
% find part of paths corresponding to relevance index
sparse_indA = [sparse_ind(sparse_ind < relevance_ind(1)),relevance_ind(1)];
sparse_indB = [sparse_ind(sparse_ind < relevance_ind(2)),relevance_ind(2)];
pathA.pos_sparse = posA(sparse_indA);
pathB.pos_sparse = posB(sparse_indB);
% sparse indeces for relevant (and sparse) paths
pathA.s_ind = sparse_indA;
pathB.s_ind = sparse_indB;
% collect into two fields that contains indeces and paths
paths.A = pathA;
paths.B = pathB;

% number of positions in the sparse paths (used mainly for plotting)
paths.n_sparse     = [length(pathA.s_ind), length(pathB.s_ind)];
% true number of positions of relevant paths (not sparse)
paths.n_rel        = relevance_ind;
% proximity time approximations
paths.t_prox       = ind2time(candidates.prime_cand);
% minimum time appriximations
paths.t_min        = ind2time(i_min);
paths.contact      = contact;
% time of first path overlap (if there was one)
paths.t_first_cont = ind2time(i_first_cont);
% amount of time corresponding to relevant part of paths
paths.rel_time     = ind2time([sparse_indA(end),sparse_indB(end)]);
% slope of paths at the end of relevant path
paths.theta       = [stateA.theta,stateB.theta];

end