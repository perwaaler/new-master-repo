function time_sep = temporal_sep(enc,RUprop,max_delta,stat,plots,...
                                 w,max_iter,time_travel, path_div)
% simulates predicted paths for each road user and returns which road user
% arrives first at the collision point along with time separation and time
% advantage.

% set default values
if nargin == 5
    w = 10;
    max_iter = 80;
    time_travel = 6;
    path_div = 20;
    predpath = 1;
    
elseif nargin == 6
    max_iter = 80;
    time_travel = 6;
    path_div = 20;
    predpath = 1;
    
elseif nargin == 7
    time_travel = 6;
    path_div = 20;
    predpath = 0;
    
elseif nargin == 8
    path_div = 20;
    predpath = 1;
end

% in case interaction starts early
if stat.frame <= time_travel
    time_travel = stat.frame - 1;
end

% define structure that contains previous time_travel positions and current
% position
prev_states = enc.states(stat.frame-time_travel:stat.frame); 
Z = prev_states{1};

% variables that save positions at each iteration
pathA = NaN(1,500);
pathB = NaN(1,500);

% structure that saves field states for each iteration
ZZ = cell(1,500); 

% initiate variables used to end loop if paths diverge
D0 = 100;
tik = 0;

% initiate
ZZ{1}    = Z;
candidatesB = 1;
candidatesA = 1;

for k=1:max_iter
    
    %%%% path A
    % remove elements with indeces that spill over 
    candidatesB = candidatesB(candidatesB>=1);
    candidatesB = candidatesB(candidatesB<=k);
    
    % find distance between last position of A and candidates of B
    diffA = pathA(k) - pathB(candidatesB);
    diffA = sqrt(diffA.*conj(diffA));
    
    % update closest position in B
    [min_distA,closest{2}] = min(diffA);
    closest{2} = candidatesB(closest{2});

    % recenter around new candidate
    candidatesB = closest{2}-w:closest{2}+w;
    
    %%%% path B
    % remove elements with indeces that spill over 
    candidatesA = candidatesA(candidatesA>=1);
    candidatesA = candidatesA(candidatesA<=k);
    
    % find distance between last position of A and candidates of B
    diffB = pathB(k) - pathA(candidatesA);
    diffB = sqrt(diffB.*conj(diffB));
    
    % update closest position in B
    [min_distB, closest{1}] = min(diffB);
    closest{1} = candidatesA(closest{1});

    % recenter around new candidate
    candidatesA = closest{1}-w:closest{1}+w;
    
    % RU2 is the RU that arrives latest at the collision point
    [D,RU2] = min([min_distA, min_distB]);
    % check if paths are approaching each other
    if D0 - D < 1e-4
        tik = tik + 1
        pause(3)
        
        if tik > path_div
            k = max_iter; %#ok<*NASGU,FXSET>
            break
        end
    end
    
    if D < 2*RUprop.r
        k
        % first road user to reach collision point
        RU1 = 3 - RU2;
        steps_taken = k - 1;
        
        % extract states for time TTPO and TTPC
        F(1) = ZZ{closest{RU1}}(RU1);
        F(2) = ZZ{k}(RU2);
        % reverse time to find exact moment of collision
        t = calc_ttc(F,RUprop.r);
        
        % Time To Potential Collision
        TTPC =  steps_taken - time_travel + t;
        % Time To Path Overlap
        TTPO = closest{RU1} + t - time_travel;
        % Time advantage
        Tadv = TTPC - TTPO;
        
        % save variables
        time_sep.RU1    = RU1;
        time_sep.TTPC   = TTPC;
        time_sep.TTPO   = TTPO;
        time_sep.Tadv   = (-1)^(RU1+1)*Tadv;
        break
    end
    
    %%%% take new step
    % take step and update states
    if k > time_travel
        
        % generate new step A
        if Z(1).EA == 0
            % take NEA step
            E_speedA = normrnd(E_speed_A(Z, RUprop), RUprop.Espeed_std(1));
            E_thetaA = normrnd(E_theta_A(Z, RUprop), RUprop.Etheta_std(1));
            Z(1) = take_step(Z(1), E_speedA, E_thetaA, max_delta(1)); 
        else
            % predict step assuming constant behaviour
            Z(1) = take_step(Z(1));
        end
        
        % generate new step B
        if Z(2).EA == 0
            % take NEA step
            E_speedB = normrnd(RUprop.avg_speed(2),  RUprop.Espeed_std(2));
            E_thetaB = normrnd(pi,                   RUprop.Etheta_std(2));
            Z(2) = take_step(Z(2), E_speedB, E_thetaB, max_delta(2));
        else
            % predict step assuming constant behaviour
            Z(2) = take_step(Z(2));
        end
    end
    
    % save positions
    pathA(k+1) = Z(1).pos;
    pathB(k+1) = Z(2).pos;

    % save state
    ZZ{k+1}    = Z;
    
    % plot prdicted paths
    if plots.predpath==1
        hold on
        plot_pos(Z, plots,plots.col)
    end
end

if D > 2*RUprop.r && k==max_iter
    time_sep.RU1  = nan;
    time_sep.TTPO = nan;
    time_sep.Tadv = nan;
    time_sep.TTPC = nan;
end



end





