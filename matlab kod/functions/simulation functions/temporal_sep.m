function time_sep = temporal_sep(enc,RUprop,max_delta,stat,plots,...
                                 w, max_iter, t_prev, path_div)
% simulates predicted paths for each road user and returns which road user
% arrives first at the collision point along with time separation and time
% advantage.

% set default values
if nargin == 5
    w = 8;
    max_iter = 80;
    t_prev = 6;
    path_div = 20;
    predpath = 1;
    
elseif nargin == 6
    max_iter = 80;
    t_prev = 6;
    path_div = 20;
    predpath = 1;
    
elseif nargin == 7
    t_prev = 6;
    path_div = 20;
    predpath = 0;
    
elseif nargin == 8
    path_div = 20;
    predpath = 1;
end

% in case interaction starts early
if stat.frame <= t_prev
    t_prev = 0;
end

% define structure that contains states from time t_curr-t_prev up to 
% and including t_curr
%%
r = RUprop.r;
prev_states = enc.states(stat.frame - t_prev:stat.frame); 

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
D0 = 100;
tik = 0;
% used to see if paths are approaching
prev_min_dist = [100,100];
first_div = [0,0];

% initiate
Z = prev_states{1};
ZZ{1}       = Z;
candidatesB = 1;
candidatesA = 1;
candidates.index = {1,1};

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

    % if RU starts getting further away from other RUs trajectory, save
    % information for classification of encounter
    if D_diff(1) >= 0 && first_div(1)==0
        divA_candidates = candidates;
        first_div(1)= k;
        D_divA = candidates.min_dist(1);
    end
    if D_diff(2) >= 0 && first_div(2)==0
        divB_candidates = candidates;
        first_div(2) = k;
        D_divB = candidates.min_dist(2);
    end
    
    % the following condition indicates paths will propbably never collide
    if min(D_diff)>=0 && k>t_prev && D>2*r 

        [~, who_moves] = max([Z(1).speed, Z(2).speed]);
        
        if who_moves == 1
            m = first_div(1);
            A_pert = ZZ{m-1}(1).pos + ZZ{m}(1).speed*exp(1i*(ZZ{m}(1).theta + deg2rad(10) ));
            pathA_pert    = pathA;
            pathA_pert(m) = A_pert;
            candidates_pert = find_closest_ball(divA_candidates,...
                                                pathA_pert,pathB,m,w);
            D_diff_pert = candidates_pert.min_dist(1) - D_divA;
        else
            m = first_div(2);
            B_pert = ZZ{m-1}(2).pos + ZZ{m}(2).speed*exp(1i*(ZZ{m}(2).theta + deg2rad(10) ));
            pathB_pert    = pathB;
            pathB_pert(m) = B_pert;
            candidates_pert = find_closest_ball(divB_candidates,...
                                                pathA,pathB_pert,m,w);
            D_diff_pert = candidates_pert.min_dist(2) - D_divB;
        end

        tik = tik + 1;
        if tik>3
            TTPC = abs(enc.states{stat.frame}(1).pos - enc.states{stat.frame}(2).pos)...
                /(enc.states{stat.frame}(1).speed + enc.states{stat.frame}(2).speed);

            if D_diff_pert < 0
                Tadv=-inf;
            else
                Tadv=inf;
            end
            % save variables
            time_sep.RU1    = nan;
            time_sep.TTPC   = TTPC;
            time_sep.T2     = inf;
            time_sep.Tadv   = Tadv;
            break
        end
    end

    % check if paths are approaching each other
    if D0 - D < 1e-4
        tik = tik + 1;
        if tik > path_div
            k = max_iter; %#ok<*NASGU,FXSET>
            break
        end
    end
    
    % if this condition is satisfied, stop generating new steps
    if D > 2*RUprop.r && first_contact==0
        overlapRU1 = find_overlap_RU1(k,pathA,pathB,RU1,overlapRU2_ind,...
                                            ind_first_touch,RUprop,3);
        Tadv = overlapRU2_ind(1)-overlapRU1(1);
        TTPC = overlapRU1(1) - 2 - t_prev;
        T2   = overlapRU2_ind(1) - 2 - t_prev;
        % save variables
        time_sep.RU1    = RU1;
        time_sep.TTPC   = TTPC;
        time_sep.T2     = T2;
        time_sep.Tadv   = (-1)^(RU1+1)*Tadv;
        
%         F(1) = ZZ{overlapRU1(1)}(RU1);
%         F(2) = ZZ{k}(RU2);
%         % reverse time to find exact moment of collision
%         t = calc_ttc(F,RUprop.r)
        
        break
    end
    
    if D < 2*RUprop.r
        % save coordinates for first contact
        if first_contact == 1
            RU1 = 3 - RU2;
            ind_first_touch = [candidates.closest{RU1},candidates.closest{RU2}];
            overlapRU2_ind = k;
            first_contact = 0;
        % if not first contact, save coordinates for overlap of RU1
        else
            overlapRU2_ind(end+1) = overlapRU2_ind(end) + 1; %#ok<*AGROW>
        end
    end
    

%     if D < 2*RUprop.r
%         
%         % first road user to reach collision point
%         RU1 = 3 - RU2;
%         steps_taken = k - 1;
% 
%         % extract states for time TTPO and TTPC
%         F(1) = ZZ{candidates.closest{RU1}}(RU1);
%         F(2) = ZZ{k}(RU2);
%         % reverse time to find exact moment of collision
%         t = calc_ttc(F,RUprop.r);
% 
% 
%         % Time To Path Overlap
%         TTPO =   steps_taken - t_prev + t - 1;
%         % Time To Potential Collision
%         TTPC = candidates.closest{RU1} + t - (t_prev+1);
%         % Time advantage
%         Tadv = TTPO - TTPC;
%         
% 
%         % save variables
%         time_sep.RU1    = RU1;
%         time_sep.TTPC   = TTPC;
%         time_sep.TTPO   = T2;
%         time_sep.Tadv   = (-1)^(RU1+1)*Tadv;
%         break
%     end

    prev_min_dist = candidates.min_dist;

    %%%% take new step
    % take step and update states

    if k > t_prev

        % generate new step A
        if Z(1).EA == 0
            % take NEA step
            E_speedA = normrnd(E_speed_A(Z, RUprop), RUprop.Espeed_std(1));
            E_thetaA = normrnd(E_theta_A(Z, RUprop), RUprop.Etheta_std(1));
            Z(1) = take_step(Z(1), E_speedA, E_thetaA, max_delta(1)); 
        else
%             predict step assuming constant behaviour
            Z(1) = take_step(Z(1));
        end

%         generate new step B
        if Z(2).EA == 0
            % take NEA step
            E_speedB = normrnd(RUprop.avg_speed(2),  RUprop.Espeed_std(2));
            E_thetaB = normrnd(pi,                   RUprop.Etheta_std(2));
            Z(2) = take_step(Z(2), E_speedB, E_thetaB, max_delta(2));
        else
%             predict step assuming constant behaviour
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

% if max iterations has been reached, and there loop has still not reached
% termination criteria
if k==max_iter
    time_sep.RU1  = nan;
    time_sep.T2   = nan;
    time_sep.Tadv = nan;
    time_sep.TTPC = nan;
end

end


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
    plot_pos(ZZ{overlap(i)}, plots)
end
pause(1)
end


%%
end