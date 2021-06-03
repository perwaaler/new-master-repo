function states = gen_random_states()

speed0 = max(0.01, gam_rnd(0.20, 0.05, 2)); 
speed  = speed0; %#ok<*NASGU>
theta  = normrnd([0,pi], deg2rad(7));
dtheta = [0,0];
dspeed = [0,0];
% determine initial x values
init_x = 6;
A_real = -6;
B_real = 6;
A_im = normrnd(-2.3,0.5);
B_im = normrnd(0,   0.5);
A0 = A_real + 1i*A_im;
B0 = B_real + 1i*B_im;
% set initial decision states
decision = [0,0];

% choose how similar a trajectory has to be in order to be "similar"
tolerance = 0.4*[0.16*0.1, 0.03, 0.08]; 


% set regression parameters for probability to engage in Evasive Action(EA)
beta_EA.int   = -2.8;
beta_EA.dist  = .8;
beta_EA.speed = -1;

%%%% set driver properties %%%%
% radius of Road Users
RUprop.r         = 0.3;
% determine how jerky the RU are in their movements
RUprop.Etheta_std      = beta_rnd(deg2rad(3),deg2rad(.5),2); 
RUprop.Espeed_std      = beta_rnd(0.05, 0.01,2);
RUprop.d2theta_std     = beta_rnd(deg2rad(7),deg2rad(2),2);
% set how accurately the drivers can predict future positions
RUprop.pred_speed_err = gam_rnd(0.001, 0.0001,2);
RUprop.pred_theta_err = gam_rnd(deg2rad(.5),deg2rad(.05),2);
% set minimum speed
RUprop.speed_min = 0.01;
RUprop.avg_speed = speed0;
% determine how much A brakes during left turn
RUprop.turn_brake = beta_rnd(0.03,0.003);
RUprop.start_brake = beta_rnd(0.03,0.003);
% weight determines momentum of drivers
RUprop.weight    = gam_rnd(1.0,0.09,2);
% where does A start and end turn
RUprop.turnstart = normrnd(-5.2,0.5);
RUprop.turnend   =  RUprop.turnstart + gam_rnd(5,0.5);
% determine how alert each driver is
RUprop.alert = beta_rnd(.8,.08,2);
% determine how often RU assesses success of decision
RUprop.decision_freq = [2,2];
RUprop.aggression = beta_rnd(0.5,0.02,2);
RUprop.id = [1,2];

%%%% set laws of movement / laws of momentum %%%%
max_delta(1).dspeed = beta_rnd(0.01*RUprop.weight(1)^-1, 1e-4);
max_delta(2).dspeed = beta_rnd(0.01*RUprop.weight(2)^-1, 1e-4);
max_delta(1).dtheta = beta_rnd(0.3*RUprop.weight(1)^-1, 8e-4);
max_delta(2).dtheta = beta_rnd(0.3*RUprop.weight(2)^-1, 8e-4);
max_delta(1).d2theta = deg2rad(15);
max_delta(2).d2theta = deg2rad(15);
% for convenience set laws of physics as a driver property
RUprop.max_delta = max_delta;


% generate 2 paths:
stateA.pos    = normrnd(-2,1) + normrnd(-2.5,2)*1i;
stateA.theta  = normrnd(deg2rad(10),deg2rad(7));
stateA.dtheta = normrnd(deg2rad(10), deg2rad(6));
stateA.speed  = max([0,gam_rnd(0.3,0.02)-0.05]);
stateA.dspeed = normrnd(0,0.03);  
stateA.RUprop = RUprop;
stateA.id     = 1;

stateB.pos    = normrnd(2,1) + normrnd(1,1)*1i;
stateB.theta  = normrnd(pi,deg2rad(10));
stateB.dtheta = normrnd(0,deg2rad(3));
stateB.speed  = max([0,gam_rnd(0.3,0.02)-0.1]);
stateB.dspeed = gam_rnd(0.01,0.001);
stateB.RUprop = RUprop;
stateB.id     = 2;
plotit        = 0;
pauseL        = 2;

states{1} = stateA;
states{2} = stateB;
end






