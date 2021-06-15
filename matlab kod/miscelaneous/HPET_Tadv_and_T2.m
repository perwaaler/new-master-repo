pathA = enc.posA(real(enc.posA)<inf);
pathB = enc.posB(real(enc.posB)<inf);
clf
plot(pathA(1:end))
hold on
plot(pathB(1:end))
ylim([-3,3])
hold off
%%

%%%% initiate encounter %%%%
% determines average and initial speed of each vehicle
speed0 = max(0.01, gam_rnd(0.26, 0.05, 2)); speed  = speed0;
theta  = normrnd([0,pi], deg2rad(7));
dtheta = [0,0];
% determine initial x values
init_x = 6;
A_real = -init_x;
B_real = init_x;
A_im = normrnd(-2.3,0.5);
B_im = normrnd(0,   0.5);
A0 = A_real + 1i*A_im;
B0 = B_real + 1i*B_im;

pathA = NaN(500,1);
pathB = NaN(500,1);
pathA(1) = A0;
pathB(1) = B0;
%%
%%%% set state of each Road User (RU) %%%%

% neighbourhood size and disc radius
w = 7;
r = 0.3;
max_iter = 30;
time_travel = 6;
path_div = 5;
% define structure that contains previous time_travel positions and current
% position
prev_states = enc.states(stat.frame-time_travel:stat.frame); 
Z = prev_states{1};

% variables that save positions at each iteration
pathA = NaN(1,500);
pathB = NaN(1,500);

% structure that saves field states for each iteration
ZZ = cell(1,500); 
D0 = 100;
tik = 0;

% initiate
ZZ{1}    = Z;
candidatesB = 1;
candidatesA = 1;
clf
plot_pos(Z, plots, init_x,["blue","green"])

for k=1:max_iter

    %%%% path B
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
    
    % RU2 is the RU that arrives lates at the collision point
    [D,RU2] = min([min_distA, min_distB]);
    
    % check if paths are approaching each other
    if D > D0 + 1e-6
        tik = tik + 1;
        if tik > path_div
            k = max_iter; %#ok<*FXSET>
            break
        end
    end
    pause(.1)
    
    if D < 2*r
        steps_taken = k-1;
        RU1 = 3-RU2;
        
        % extract states from list of saved states
        F(1) = ZZ{closest{RU1}}(RU1);
        F(2) = ZZ{k}(RU2);
        t = calc_ttc(F,RUprop.r);

        
        % first road user to enter collision point
        RU1 = 3-RU2;
        time_separation = steps_taken - time_travel + t;
        % time to collision pointTTCP = closest{RU1} + t - time_travel;
        % Time advantage
        Tadv = time_separation - TTCP;
        break

    end
    
    D0 = D;
    
%     %%%% take new step
%     E_speedA = normrnd(E_speed_A(Z, RUprop), RUprop.speed_std(1));
%     E_thetaA = normrnd(E_theta_A(Z, RUprop), RUprop.theta_std(1));
%     E_speedB = normrnd(RUprop.avg_speed(2),  RUprop.speed_std(2));
%     E_thetaB = normrnd(pi,                   RUprop.theta_std(2));
%     
%     % take step and update states
%     if k > time_travel
%         % generate new step
%         Z(1) = take_step(Z(1), E_speedA, E_thetaA, max_delta(1)); 
%         Z(2) = take_step(Z(2), E_speedB, E_thetaB, max_delta(2));
%     else
%         Z(1) = prev_states{k+1}(1);
%         Z(2) = prev_states{k+1}(2);
%     end
    
        % take step and update states
    if k > time_travel
        
        % generate new step A
        if Z(1).EA == 0
            % take NEA step
            E_speedA = normrnd(E_speed_A(Z, RUprop), RUprop.speed_std(1));
            E_thetaA = normrnd(E_theta_A(Z, RUprop), RUprop.theta_std(1));
            Z(1) = take_step(Z(1), E_speedA, E_thetaA, max_delta(1)); 
        else
            % predict step assuming constant behaviour
            Z(1) = take_step(Z(1));
        end
        
        % generate new step B
        if Z(2).EA == 0
            % take NEA step
            E_speedB = normrnd(RUprop.avg_speed(2),  RUprop.speed_std(2));
            E_thetaB = normrnd(pi,                   RUprop.theta_std(2));
            Z(2) = take_step(Z(2), E_speedB, E_thetaB, max_delta(2));
        else
            % predict step assuming constant behaviour
            Z(2) = take_step(Z(2));
            Z(2).pos;
        end
    else
        Z(1) = prev_states{k+1}(1);
        Z(2) = prev_states{k+1}(2);

    end
    
    % save positions
    pathA(k+1) = Z(1).pos;
    pathB(k+1) = Z(2).pos;

    hold on
    plot_pos(Z, plots, init_x,["blue","green"])

    % save state
    ZZ{k+1}    = Z;

    
end
%%
hold on
plot_pos(Z,D, plots, init_x, r,["black","black"])

%%

%
clf
m = k;
if RU1 == 1 % A is first to arrive
    plot(pathA(1:m))
    hold on
    plot(r*exp(1i*linspace(0,2*pi,100)) + pathA(TTCP),'green')
    plot(r*exp(1i*linspace(0,2*pi,100)) + pathA(k),'red')
    
    plot(pathB(1:m))
    plot(r*exp(1i*linspace(0,2*pi,100)) + pathB(TTCP),'green')
    plot(r*exp(1i*linspace(0,2*pi,100)) + pathB(k),'red')
    ylim([-4,4])

else        % B is first to arrive
    
    plot(pathA(1:m))
    hold on
    plot(r*exp(1i*linspace(0,2*pi,100)) + pathA(TTCP),'green')
    plot(r*exp(1i*linspace(0,2*pi,100)) + pathA(k-1),'red')
    plot(r*exp(1i*linspace(0,2*pi,100)) + pathA(k),'red')
    
    plot(pathB(1:m))
    plot(r*exp(1i*linspace(0,2*pi,100)) + pathB(TTCP-1),'green')
    plot(r*exp(1i*linspace(0,2*pi,100)) + pathB(TTCP),'green')
    plot(r*exp(1i*linspace(0,2*pi,100)) + pathB(k),'red')
    ylim([-4,4])
    hold off

end
