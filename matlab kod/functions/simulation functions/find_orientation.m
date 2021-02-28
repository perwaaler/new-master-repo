function side = find_orientation(RUinfo,ZZ,RUprop)
% find orientation of RU's. RUinfo contains information about the who is
% RU1 and the time-id of the balls. ZZ is a cell that contains the
% dataframes for the whole encounter.
% RUinfo = [RU1, RU1_pos, RU2_pos].

max_delta = RUprop.max_delta;
% extract information
RU1   = RUinfo(1);
idRU1 = RUinfo(2);
idRU2 = RUinfo(3);

RU2 = 3-RU1;

states{RU2} = ZZ{idRU2}(RU2);
states{RU1} = ZZ{idRU1}(RU1);

stateA = states{1};
stateB = states{2};

if stateA.speed==0
    E_speedA = normrnd(E_speed_A(stateA, RUprop), RUprop.Espeed_std(1));
    E_thetaA = normrnd(E_theta_A(stateA, RUprop), RUprop.Etheta_std(1));
    stateAdesired = take_step(stateA, E_speedA, E_thetaA, max_delta(1));
    A0 = stateA.pos;
    A1 = stateAdesired.pos;
else
    A0 = stateA.pos - stateA.speed*exp(1i*stateA.theta);
    A1 = stateA.pos;
end

if stateB.speed==0
    E_speedB = normrnd(RUprop.avg_speed(2),  RUprop.Espeed_std(2));
    E_thetaB = normrnd(pi,                   RUprop.Etheta_std(2));
    % take step and update states
    stateBdesired = take_step(stateB, E_speedB, E_thetaB, max_delta(2));
    B0 = stateB.pos;
    B1 = stateBdesired.pos;
else
    B0 = stateB.pos - stateB.speed*exp(1i*stateB.theta);
    B1 = stateB.pos;
end
    
% A0 = stateA.pos - stateA.speed*exp(1i*stateA.theta);
% A1 = stateA.pos;
% B0 = stateB.pos - stateB.speed*exp(1i*stateB.theta);
% B1 = stateB.pos;

% define tangent vectors
if RU1==1
    u = A1-A0;
    v = B1-A0;
else
    u = B1-B0;
    v = A1-B0;
end

side = find_side(u,v);

% 
% plot_circle(A0,r,"magenta")
% plot_circle(B1,r,"cyan")
% plot_circle(A1,r,"black")


end
