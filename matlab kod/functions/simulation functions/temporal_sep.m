function time_sep = temporal_sep(enc,RUprop,stat,plots,...
                                 w, max_iter, t_prev, divtol)
% simulates predicted paths for each road user and returns which road user
% arrives first at the collision point along with time separation and time
% advantage.

% extract variables from fields
r = RUprop.r;
max_delta = RUprop.max_delta;

% set default values
if nargin == 4
    w = 8;
    max_iter = 80;
    t_prev = 6;
    divtol = 3;
elseif nargin == 5
    max_iter = 80;
    t_prev = 6;
    divtol = 3;
elseif nargin == 6
    t_prev = 6;
    divtol = 3;
elseif nargin == 7
    divtol = 3;
end

% in case interaction starts early
if stat.frame <= t_prev
    t_prev = 0;
end


%%
% clf
% stat.frame=stat.frame-1
% plots.predpath=1
% time_sep

% define structure that contains states from time t_curr-t_prev up to 
% and including t_curr
prev_states = enc.states(stat.frame - t_prev:stat.frame); 
state0 = prev_states{end};

% average out dtheta (turn rate) to make extrapolations more stable
if t_prev > 0
    dtheta = [0,0];
    for i = 1:t_prev
        dtheta(1) = dtheta(1) + prev_states{i}(1).dtheta;
        dtheta(2) = dtheta(2) + prev_states{i}(2).dtheta;
    end
    dtheta = dtheta/t_prev;
    prev_states{stat.frame}(1).dtheta = dtheta(1);
    prev_states{stat.frame}(2).dtheta = dtheta(2);
end

% variables that save positions at each iteration
pathA = NaN(1,500);
pathB = NaN(1,500);

pathA(1) = prev_states{1}(1).pos;
pathB(1) = prev_states{1}(2).pos;

% structure that saves field states for each iteration
ZZ = cell(1,500); 

% initiate variables used to end loop if paths diverge
tik = 0;
% used to see if paths are approaching
prev_min_dist = [100,100];

% initiate
Z = prev_states{1};
ZZ{1}       = Z;
candidates.index = {1,1};
speeds = [ZZ{1}(1).speed,ZZ{1}(2).speed];
[~, id_fastest] = max(speeds);

% indicates whether or not path-overlap has occured
first_contact = 1;
% saves index of balls that overlap with atleast one ball in path of RU1
overlapRU2_ind = [];

for k=1:max_iter

    candidates = find_closest_ball(candidates,...
                                            pathA, pathB, k, w);
    % RU2 is the RU that arrives latest to the collision point
    [D,RU2] = min(candidates.min_dist);
    
    % are paths getting closer? Positive D_diff indicates growing distance
    D_diff = candidates.min_dist - prev_min_dist;
    % set indicator for divergence status
    paths_are_div = min(D_diff)>=0;
    prev_min_dist = candidates.min_dist;

   
    % used to end loop if paths are not moving
    speeds = [Z(1).speed,Z(2).speed];
    if max(speeds)==0
        % Both have stopped. Let id_fastest be whatever it was in 
        % previous iteration. Set tik=3 to make sure loop gets ended
        tik = divtol;
    else
        [~, id_fastest] = max(speeds);
        one_RU_has_stopped = min(speeds)==0;
    end
    
    
    % EXIT 1: PATH OVERLAP
    % this condition implies that the overlap of RU2 has be computed, and 
    % that the overlap of RU1 can be calculated
    if (D > 2*r || k==max_iter ) && first_contact==0
        overlapRU1 = find_overlap_RU1(k,pathA,pathB,RU1,overlapRU2_ind,...
                                            ind_ft,RUprop,3);
        Tadv = overlapRU2_ind(1)-overlapRU1(1);
        TTPC = overlapRU1(1) - 2 - t_prev;
        T2   = overlapRU2_ind(1) - 2 - t_prev;
        
        % save variables
        time_sep.RU1    = RU1;
        time_sep.TTPC   = TTPC;
        time_sep.T2     = T2;
        time_sep.Tadv   = (-1)^(RU1+1)*Tadv;
        time_sep.side   = [side,side];
        break
    end
    
    
    % EXIT 2: NO PATH OVERLAP
    % the following event indicates paths are diverging AND not overlapping
    if (paths_are_div==1  || k==max_iter) && D>2*r && k>t_prev
        % whoever moves fastest is determined to be RU1
        RU1 = id_fastest;  
        tik = tik + 1;
        
        if tik > divtol
            % set TTPC to be TTC from current moment if they were heading
            % straight towards eachother with unchanged speed.
            TTPC = abs(enc.states{stat.frame}(1).pos - enc.states{stat.frame}(2).pos)...
                /(enc.states{stat.frame}(1).speed + enc.states{stat.frame}(2).speed);
            
            if one_RU_has_stopped==1
                n_samplePts = 7;
            else
                n_samplePts = 5;
            end
            
            Tadv = inf;
            side = find_avg_orientation(ZZ,RU1,k,t_prev,n_samplePts,RUprop);

            % save variables
            time_sep.RU1  = RU1;
            time_sep.TTPC = TTPC;
            time_sep.T2   = inf;
            time_sep.Tadv = Tadv;
            time_sep.side = side;
            break
        end
    end
    
    
    % computes the overlap of RU2 once overlapp has begun
    if D < 2*r
        % save coordinates of the moment of first contact
        if first_contact == 1
            RU1 = 3 - RU2;
            % If RU1 ambiguous, use linear extrapolation to find RU1:
            if candidates.closest{1}==candidates.closest{2}
                RU1 = RU1_tiebraker(Z,r);
                RU2 = 3-RU1;
            end
            
            % save index for each ball at first touch
            ind_ft = candidates.neighbours(RU2,:);
            ind_ft = [ind_ft(2),ind_ft(1)];
            % collect first overlap index of RU2
            overlapRU2_ind = k;
            first_contact = 0;
            % find which side of the of the path of RU1 that RU2 has hit
            side = find_orientation([RU1,ind_ft],ZZ,RUprop);
            
        % if not first contact, save coordinates for overlap of RU1
        else
            overlapRU2_ind(end+1) = overlapRU2_ind(end) + 1; %#ok<*AGROW>
        end
    end

    
    %%%% take new step %%%%
    % take step and update states:
    if k > t_prev

        % generate new step A
        if Z(1).EA == 0
            % take NEA step
            E_speedA = normrnd(E_speed_A(Z, RUprop), RUprop.Espeed_std(1));
            E_thetaA = normrnd(E_theta_A(Z, RUprop), RUprop.Etheta_std(1));
            Z(1) = take_step(Z(1), E_speedA, E_thetaA, max_delta(1)); 
        else
            % straighten out to prevent going in circles
            if abs(Z(1).theta - state0(1).theta) > pi
                Z(1).dtheta = 0;
            end
%           predict step assuming constant behaviour
            Z(1) = take_step(Z(1));
        end

%         generate new step B
        if Z(2).EA == 0
            % take NEA step
            E_speedB = normrnd(RUprop.avg_speed(2),  RUprop.Espeed_std(2));
            E_thetaB = normrnd(pi,                   RUprop.Etheta_std(2));
            Z(2) = take_step(Z(2), E_speedB, E_thetaB, max_delta(2));
        else
            % straighten out to prevent going in circles
            if abs(Z(2).theta - state0(2).theta) > pi
                Z(2).dtheta = 0;
            end
%           predict step assuming constant behaviour
            Z(2) = take_step(Z(2));
        end
        
    else
        Z = prev_states{k+1};
    end

    
    % save positions
    pathA(k+1) = Z(1).pos;
    pathB(k+1) = Z(2).pos;

    % save state
    ZZ{k+1}    = Z;

    
    % plot predicted paths
    if plots.predpath==1
        if k==t_prev+1
            plots.col = ["black","black"];
            hold on
            plot_pos(Z, plots)
        else
            plots.col = ["blue","red"];
            hold on
            plot_pos(Z, plots)    
        end        
    end

end


% in case we want to do some plotting of the predicted paths
if plots.overlap==1 && plots.predpath==1 && first_contact==0 
common_val = intersect(overlapRU1,overlapRU2_ind);
overlap = [overlapRU1,overlapRU2_ind];
% determines color-scheme for the overlapping balls
col_matrix = [["blue","green"];
              ["red","yellow"]];

for i=1:length(overlap)
    
    if     ismember(overlap(i),common_val)
        plots.col = col_matrix(:,2);
        
    elseif ismember(overlap(i),overlapRU1)
        plots.col = [col_matrix(1,RU2),col_matrix(2,RU1)];
        
    elseif ismember(overlap(i),overlapRU2_ind)
        plots.col = [col_matrix(1,RU1),col_matrix(2,RU2)];
    end
    
    hold on
    plot_pos(ZZ{overlap(i)}, plots, time_sep)
end
pause(1)
end

%%
end















% OLD CODE %%%%%%%%%%%%%

% USED TO IDENTIFY MOMENT OF CLOSEST PROXIMITY
%
%     % this information is used to identify moment of least separation as
%     % well as the index of the road users at that frame. Used to classify
%     % orientation when there is no path overlap.
%     if D < track_closest
%         % in the no-overlap-case, RU1 is defined to be the fastest RU
%         RU1no = id_fastest;
%         RU2no = 3 - RU1no;
%         closest_neighbours = candidates.neighbours(RU1no,:);
%         save_closest = [RU1no, closest_neighbours];
%         track_closest = D;
%     end

% CODE USED TO CHARACTERIZE ORIENTATION USING ANGLE PERTURBATIONS
%
%     % if RU starts getting further away from other RUs trajectory, save
%     % information for classification of encounter
%     if D_diff(1) >= 0 && t_first_div(1)==0
%         divA_candidates = candidates;
%         % identifies moment of first divergence
%         t_first_div(1)= k;
%         D_divA = candidates.min_dist(1);
%     end
%     if D_diff(2) >= 0 && t_first_div(2)==0
%         divB_candidates = candidates;
%         t_first_div(2) = k;
%         D_divB = candidates.min_dist(2);
%     end

%         if id_fastest == 1
%             m = t_first_div(1);
%             A_pert = ZZ{m-1}(1).pos + ZZ{m}(1).speed*exp(1i*(ZZ{m}(1).theta + deg2rad(10) ));
%             pathA_pert    = pathA;
%             pathA_pert(m) = A_pert;
%             candidates_pert = find_closest_ball(divA_candidates,...
%                                                 pathA_pert,pathB,m,w);
%             D_diff_pert = candidates_pert.min_dist(1) - D_divA;
%         else
%             m = t_first_div(2);
%             B_pert = ZZ{m-1}(2).pos + ZZ{m}(2).speed*exp(1i*(ZZ{m}(2).theta + deg2rad(10) ));
%             pathB_pert    = pathB;
%             pathB_pert(m) = B_pert;
%             candidates_pert = find_closest_ball(divB_candidates,...
%                                                 pathA,pathB_pert,m,w);
%             D_diff_pert = candidates_pert.min_dist(2) - D_divB;
%         end
