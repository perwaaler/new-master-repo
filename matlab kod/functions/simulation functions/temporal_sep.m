function time_sep = temporal_sep(S,RUprop,max_delta,search_width,max_iter)
% simulates predicted paths for each road user and returns which road user
% arrives first at the collision point along with time separation and time
% advantage.

% set default values
if nargin == 3
    search_width = 3;
    max_iter = 50;
elseif nargin == 4
    max_iter = 50;
end

% variables that save positions at each iteration
pathA = NaN(1,500);
pathB = NaN(1,500);

% structure that saves field states for each iteration
SS = cell(1,500); 

% initiate
SS{1}    = S;
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
    candidatesB = closest{2}-search_width:closest{2}+search_width;
    
    %%%% path B
    % remove elements with indeces that spill over 
    candidatesA = candidatesA(candidatesA>=1);
    candidatesA = candidatesA(candidatesA<=k);
    
    % find distance between last position of A and candidates of B
    diffB = pathB(k) - pathA(candidatesA);
    diffB = sqrt(diffB.*conj(diffB));
    
    % update closest position in B
    [min_distB,closest{1}] = min(diffB);
    closest{1} = candidatesA(closest{1});

    % recenter around new candidate
    candidatesA = closest{1}-search_width:closest{1}+search_width;
    
    % RU2 is the RU that arrives latest at the collision point
    [D,RU2] = min([min_distA, min_distB]);
    
    if D < 2*RUprop.r
        % first road user to reach collision point
        RU1 = 3-RU2;
        
        % extract states for time TTCP and Tsep
        F(1) = SS{closest{RU1}}(RU1);
        F(2) = SS{k}(RU2);
        % reverse time to find exact moment of collision
        t = calc_ttc(F,RUprop.r);
        
        % Time separation
        Tsep = k + t;
        % Time To Collision Point
        TTCP = closest{RU1} + t;
        % time to collision point
        Tadv = Tsep - TTCP;
        
        % save variables
        time_sep.RU1    = RU1;
        time_sep.total  = Tsep;
        time_sep.TTCP   = TTCP;
        time_sep.Tadv   = Tadv;
        break
    end
    
    %%%% take new step
    E_speedA = normrnd(E_speed_A(S, RUprop), RUprop.speed_std(1));
    E_thetaA = normrnd(E_theta_A(S, RUprop), RUprop.theta_std(1));
    E_speedB = normrnd(RUprop.avg_speed(2),  RUprop.speed_std(2));
    E_thetaB = normrnd(pi,                   RUprop.theta_std(2));

    % take step and update states
    S(1) = take_step(S(1), E_speedA, E_thetaA, max_delta(1)); 
    S(2) = take_step(S(2), E_speedB, E_thetaB, max_delta(2));
    
    % save positions
    pathA(k+1) = S(1).pos;
    pathB(k+1) = S(2).pos;
    
    % save state
    SS{k+1}    = S;

end

if k == max_iter
    time_sep.RU1  = nan;
    time_sep.TTCP = nan;
    time_sep.Tadv = nan;
    time_sep.T2   = nan;
end

end





