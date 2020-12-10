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
[S(1).pos,    S(2).pos]       = deal(A0,B0);
[S(1).theta,  S(2).theta]     = deal(theta(1),theta(2));
[S(1).speed,  S(2).speed]     = deal(speed0(1), speed0(2));
[S(1).dtheta, S(2).dtheta]    = deal(dtheta(1),dtheta(2));
[S(1).EA,     S(2).EA]        = deal(0,0);

% neighbourhood size
w = 3;
r = 0.3;
pathA = NaN(1,500);
pathB = NaN(1,500);
SS = cell(1,500); 
SS{1}    = S;

pl = 1;
candidatesB = 1;
candidatesA = 1;
for k=1:50
    
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
    [min_distB,closest{1}] = min(diffB);
    closest{1} = candidatesA(closest{1});

    % recenter around new candidate
    candidatesA = closest{1}-w:closest{1}+w;
    
    % RU2 is the RU that arrives lates at the collision point
    [M,RU2] = min([min_distA, min_distB]);
    
    if M<2*r

        RU1 = 3-RU2;
        
        % extract states from list of saved states
        F(1) = SS{closest{RU1}}(RU1);
        F(2) = SS{k}(RU2);
        t = calc_ttc(F,RUprop.r)

        
        % first road user to enter collision point
        RU1 = 3-RU2;
        time_separation = k-t;
        % time to collision point
        TTCP = closest{RU1}
        % Time advantage
        Tadv = time_separation - TTCP
        break

    end
    
    %%%% take new step
    E_speedA = normrnd(E_speed_A(real(S(1).pos) + .3, RUprop), RUprop.speed_std(1));
    E_thetaA = normrnd(E_theta_A(S,RUprop),              RUprop.theta_std(1));
    E_speedB = normrnd(speed0(2),                        RUprop.speed_std(2));
    E_thetaB = normrnd(pi,                               RUprop.theta_std(2));

    % take step and update states
    S(1) = take_step(S(1), E_speedA, E_thetaA, max_delta(1)); 
    S(2) = take_step(S(2), E_speedB, E_thetaB, max_delta(2));
    
    pathA(k+1) = S(1).pos;
    pathB(k+1) = S(2).pos;
    
    % save state
    SS{k+1}    = S;
    
end
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
%     hold off
else
    
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
%     hold off

end
%%
    A0 = H(1,2);
    B0 = H(2,2);
    velocity = calc_speed_angle(H)
    speed0   = [velocity(1).speed,velocity(2).speed];
    theta0   = [velocity(1).theta,velocity(2).theta];
    
    H(1,1) + speed0(1)*exp(1i*theta0(1))
    H(2,1) + speed0(2)*exp(1i*theta0(2))
