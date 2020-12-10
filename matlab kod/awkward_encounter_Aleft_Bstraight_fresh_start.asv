n_est_sets = 1;                        % number of encounter-samples desired
N_enc = 100;                           % number of encounters
all_data = cell(n_est_sets, 12);       % collects data from each encounter-sample
NTTC = 30;                             % Number of TTC to sample at first evasive action for each encounter
est_ttc = 0;                           % tells the algorithm whether or not you want to estimate ttc distribution for each encounter
compute_X = 0;
[plots.enc, plots.predpath, plots.pause] =  deal(1, 0, 2^-5); % plot settings

disable_crash = 0;                       % disable collision and EA mode to record free movement patterns
use_history = 0;
use_dp_model = 0;                        % set to 1 if you want to use estimated reaction model

for kk = 1:n_est_sets
kk %#ok<NOPTS>

% EA behaviour modification parameters

% variables that collects severity meausures and data for each encounter
fea.dist = inf*ones(N_enc,1);              % saves separation distance at time of first evasive action
fea.ttc = inf*ones(N_enc,1);               % saves TTC at FEA
fea.stoch_ttc = nan*ones(N_enc, NTTC);     % saves estimates of TTC distributions; each row contains NTTC simulated ttc values

ea.mindist = inf*ones(N_enc,1000);         % saves min sep. dist. over EA part
ea.minttc = inf*ones(N_enc,1000);          % saves min TTC over EA part

enc.mindist = inf*ones(1,N_enc);           % saves min sep. dist.
enc.minttc = inf*ones(1,N_enc);            % saves min TTC
% + --> no EA, - --> EA); 1 --> no collision, 2 --> collision



for i=1:N_enc

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

% choose how similar a trajectory has to be in order to be "similar"
tolerance = 0.4*[0.16*0.1, 0.03, 0.08]; 


% set parameters for probability to engage in Evasive Action(EA)
beta_EA.int   = -2.8;
beta_EA.dist  = .8;
beta_EA.speed = -.4;

%%%% set driver properties %%%%
% radius of Road Users
RUprop.r         = 0.3;
% determine steadiness/smoothness of Road Users (RU)
RUprop.theta_std = beta_rnd(deg2rad(3),deg2rad(.5),2); 
RUprop.speed_std = beta_rnd(0.05, 0.01,2);
RUprop.speed_min = 0.01;
RUprop.avg_speed = speed0;
% determine how much A brakes during left turn
RUprop.turn_brake = beta_rnd(0.03,0.003);
RUprop.start_brake = beta_rnd(0.03,0.003);
% weight determines momentum of drivers
RUprop.weight    = gam_rnd(1,0.09,2);
% where does A start and end turn
RUprop.turnstart = normrnd(-5.2,0.5);
RUprop.turnend   =  RUprop.turnstart + gam_rnd(5,0.5);
% determine how alert each driver is
RUprop.alert = beta_rnd(.8,.08,2);                     


%%%% set laws of movement / laws of momentum %%%%
max_delta(1).dspeed = beta_rnd(0.01*RUprop.weight(1), 1e-4);
max_delta(2).dspeed = beta_rnd(0.01*RUprop.weight(2), 1e-4);
max_delta(1).dtheta = beta_rnd(0.03*RUprop.weight(1), 8e-4);
max_delta(2).dtheta = beta_rnd(0.03*RUprop.weight(2), 8e-4);
max_delta(1).d2theta = deg2rad(5);
max_delta(2).d2theta = deg2rad(5);


%%%% set state of each Road User (RU) %%%%
[S(1).pos,    S(2).pos]       = deal(A0,B0);
[S(1).theta,  S(2).theta]     = deal(theta(1),theta(2));
[S(1).speed,  S(2).speed]     = deal(speed0(1), speed0(2));
[S(1).dtheta, S(2).dtheta]    = deal(dtheta(1),dtheta(2));
[S(1).EA,     S(2).EA]        = deal(0,0);

% save state
enc.states = cell(N_enc,500); 

%%%% variables that track status %%%%
% track frame
stat.frame     = 0;
stat.EAframe   = 1;          
% have they interacted yet
stat.first_int = 0;
% is atleast one of them interacting
stat.int       = 0;

                          
D_enc_i      = inf*ones(1,500);           
D_enc_i(1)   = norm(A0-B0)-2*RUprop.r;
ttc_enc_i    = inf*ones(1,500);           % saves danger index associated with each timestep
ttc_enc_i(1) = calc_ttc(S,RUprop.r);

% if disable_crash==1
%   engage_EA = double(imag(A0)<4 || real(B0)>-xinit); 
% else
%   engage_EA = double(real(A0) < xinit   && imag(A0)<4 && real(B0) > -xinit && imag(A0)<imag(B0)+1.5 && real(A0)<real(B0)+1.5);
% end
%     while  imag(A0)<4 || real(B0)>-xinit 

    while continue_simulation(S,init_x) == 1
        
        stat.frame = stat.frame + 1;
        
        % compute probability of A or B reacting
        pea = p_interact(S,beta_EA);
  
        if S(1).EA==0 || S(2).EA==0 % if one of them has not engaged EA

            pea = p_interact(S,beta_EA);
            
            % update EA status and interaction status
            S(1).EA = max(double(rand(1)<pea(1)*RUprop.alert(1)),S(1).EA);
            S(2).EA = max(double(rand(1)<pea(2)*RUprop.alert(2)),S(2).EA);
            stat.int = max(S(1).EA,S(2).EA);
        end
        
        if disable_crash == 1 
            pea=0; 
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% DETECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % if  "drivers are in interaction mode" OR "interaction is engaged by pea roll"
        % and  "they have not yet passed each other"
        if  engage_EA(S,stat) == 1 
            
            % estimate stochastic TTC at First Evasive Action
            if stat.first_int == 0 && est_ttc == 1 
                
                if use_history==1
                    [row, col] = find_sim_sit(A0, speed(1),theta(1), historyA, tolerance);
                else
                    row=0;
                end
                
                for k=1:NTTC 

                    if k>length(row) || use_history==0 %number of observed walks satisfying init. cond.
                        ttc = ttc_simulator_double_momentumV2(S, r, driver_prop, plot_pred_path);
                        fea.stoch_ttc(i,k) = ttc;
                    else
                        ttc = ttc_simulator_history(A0,B0,speed,theta,r,driver_prop,historyA{1},row,col,k,plot_pred_path);
                        fea.stoch_ttc(i,k) = ttc;
                    end
                end
                fea.ttc(i) = calc_ttc(S,RUprop.r);
                fea.dist(i) = D;
                ea.mindist(i,stat.EAframe) = D;
                ea.minttc(i,stat.EAframe) = fea.ttc(i);
            end
            
            % update status of encounter
            stat.EAframe = stat.EAframe + 1;
            stat.first_int = 1; % set to 1 so that above loop only runs on first interaction
            S(1).EA = 1;
            enc.type(i) = -1;

            %%%% take evasive step
            E_speedA = normrnd(E_speed_A(S, RUprop), RUprop.speed_std(1));
            E_thetaA = normrnd(E_theta_A(S, RUprop), RUprop.theta_std(1));
            E_speedB = normrnd(RUprop.avg_speed(2),  RUprop.speed_std(2));
            E_thetaB = normrnd(pi,                   RUprop.theta_std(2));
            
            time_diff = temporal_sep(S,RUprop,max_delta);
            % positive Tadv means that A has a time advantage
            Tadv      = (time_diff.Tadv)*-1^(1+time_diff.RU1);
            
            if S(1).EA==1 %&& S(2).EA==0  
                p_strat = logistic_fcn()
                
                if time_diff.RU1 == 1
                    S(1) = take_step(S(1), E_speedA, E_thetaA+.5, max_delta(1));
                    S(2) = take_step(S(2), E_speedB, E_thetaB, max_delta(2));
                else
                    S(1) = take_step(S(1), 0*E_speedA, -.4, max_delta(1));
                    S(2) = take_step(S(2), E_speedB, E_thetaB, max_delta(2));
%                     pause(.1)
                    time_diff.RU1
                    time_diff.Tadv
                end
            end
            
            A0 = S(1).pos;
            B0 = S(2).pos;
            %%%%

            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% NO DETECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            % saving information in vectors, and updating status
            % find desired angle and speed of each driver
            % (add small value to argument to make driver A slow down before
            % he starts to make turn)
            E_speedA = normrnd(E_speed_A(S, RUprop), RUprop.speed_std(1));
            E_thetaA = normrnd(E_theta_A(S, RUprop), RUprop.theta_std(1));
            E_speedB = normrnd(RUprop.avg_speed(2),  RUprop.speed_std(2));
            E_thetaB = normrnd(pi,                   RUprop.theta_std(2));
            
            % take step and update states
            S(1) = take_step(S(1), E_speedA, E_thetaA, max_delta(1)); 
            S(2) = take_step(S(2), E_speedB, E_thetaB, max_delta(2));
%             pause(.5)
            A0 = S(1).pos;
            B0 = S(2).pos;
            
        end
        
        % save data from frame i of encounter kk
        D = norm(A0-B0) - 2*RUprop.r;
        D_enc_i(stat.frame) = D;
        ttc_enc_i(stat.frame) = calc_ttc(S,RUprop.r);
        
        % save state
        enc.states{i} = S;
        % plot drivers
        plot_pos(S,D, plots, init_x, RUprop.r)
        
        if D<0 % collision has occurred!
            if enc.type==-1
                % set to crash type EA
                enc.type(i) = -2;
            else
                % set to crash type no EA
                enc.type(i) = +2;
            end
            % end encounter
            break
        end
        
    end
    
    % save data from encounter kk
    enc.speedA(i) = S(1).speed;
    enc.speedB(i) = S(2).speed;

    enc.mindist(i) = min(D_enc_i);
    enc.minttc(i) = min(ttc_enc_i);

end

if disable_crash == 1
    historyA{1} = enc.posA;
    historyA{2} = enc.speedA;
    historyA{3} = enc.thetaA;
end


%plot_encounter(enc.posA, enc.posB,enc.int_stateA, enc.int_stateB, find(enc.type==-2), .15, xinit, RUprop.r)






%%%% process data

sum(enc.type==-2)
sum(enc.type==2)
% remove all elements/rows corresponding to encounters with no EA
fea.dist(enc.type==1,:) =[];
X = fea.stoch_ttc;
fea.stoch_ttc(enc.type==1,:) = [];
ea.mindist(enc.type>0,:) = [];
ea.mindist = min(ea.mindist,[],2);
ea.minttc( enc.type > 0 ,:) = [];
ea.minttc = min(ea.minttc,[],2);
dist_min_nea = dist_min_nea( enc.type > 0 );


% save data from estimation set kk


avg_cond_ttc = mean(fea.stoch_ttc,2,'omitnan');
p_hypo_coll = sum(1-isnan(fea.stoch_ttc'))'/NTTC;
weighted_ttc_fea = avg_cond_ttc.*(1+exp(-3*p_hypo_coll));
% ttc and stochastic-ttc measurements
all_data{kk,1} = fea.stoch_ttc;
all_data{kk,2} = weighted_ttc_fea;
all_data{kk,3} = fea.ttc;
all_data{kk,4} = enc.minttc';
all_data{kk,5} = ea.minttc;

% minim D measurements
all_data{kk,6} = fea.dist;
all_data{kk,7} = ea.mindist;
all_data{kk,8} = enc.mindist';

% other information
all_data{kk,9} = sum(enc.type==2);                                               % saves p_nea_coll
all_data{kk,10} = sum(enc.type==-2);                                             % saves p_ea_coll
all_data{kk,11} = (sum(enc.type==-1) + sum(enc.type==-2) + sum(enc.type==2))/N_enc;  % saves p_interactive
all_data{kk,12} = (sum(enc.type==-1)+sum(enc.type==-2))/N_enc;                       % saves p_ea


end