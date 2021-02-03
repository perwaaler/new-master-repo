function change = coulumbs_desire(S)

reactivityA = 10;
reactivityB = 10;
A0 = S(1).pos;
B0 = S(2).pos;
speedA = S(1).speed;
speedB = S(2).speed;
theta    = [S(1).theta, S(2).theta];

Aold = A0 - speedA*exp(1i*theta(1));
Bold = B0 - speedB*exp(1i*theta(2));
velA = A0-Aold;
velB = B0-Bold;
predA = A0 + [1 2]*velA;
predB = B0 + [1 2]*velB;
rA = abs(predA(1) - predB);
rB = abs(predB(1) - predA);

thetaA = angle(velA);
thetaB = angle(velB);

% calculate charge of vehicles
QA = abs(velA);
QB = abs(velB);

% calculate the force acting on each car
FA = reactivityA*QA*QB./(rA-0.2).^2 .* (predA(1)-predB);
FB = reactivityB*QA*QB./(rB-0.2).^2 .* (predB(1)-predA);

% calucate desired future positions
Adesired = predA(1) + sum(FA);
Bdesired = predB(1) + sum(FB);

% difference between where they are and where they want to be
Des_moveA = Adesired - A0;
Des_moveB = Bdesired - B0;

% calculate angle of desired move
delta_thetaA = angle(Des_moveA) - angle(predA(1)-A0);
delta_thetaB = angle(Des_moveB) - angle(predB(1)-B0);

%delta_s_max = 0.01;
delta_theta_max = 2/180*pi;

% change speed
delta_speedA = norm(Des_moveA)-norm(predA(1)-A0);
delta_speedB = norm(Des_moveB)-norm(predB(1)-B0);
% desired speed modification 
speedA_mod = sign(delta_speedA)*min(0.02, abs(delta_speedA));
speedB_mod = sign(delta_speedB)*min(0.02, abs(delta_speedB));

% mean and standard deviation of delta theta and delta speed
mu_thetaA    = delta_theta_max*(1 - exp(-1*abs(delta_thetaA)));
mu_thetaB    = delta_theta_max*(1 - exp(-1*abs(delta_thetaB)));
sigma_thetaA = deg2rad(1) + exp(-abs(delta_thetaA)*30);
sigma_thetaB = deg2rad(1) + exp(-abs(delta_thetaB)*30);
logn_parA    = lognormal_par(mu_thetaA, sigma_thetaA);
logn_parB    = lognormal_par(mu_thetaB, sigma_thetaB);

% desired angle
thetaA_mod = sign(delta_thetaA)*lognrnd(logn_parA(1), logn_parA(2));
thetaB_mod = sign(delta_thetaA)*lognrnd(logn_parB(1), logn_parB(2));

change.speed = [speedA_mod, speedB_mod];
change.theta = [thetaA_mod, thetaB_mod];

end