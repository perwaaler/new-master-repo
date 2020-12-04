n = 1;                                 % number of encounter-samples desired
N = 2000;                               % number of encounters
all_data = cell(n, 12);                % collects data from each encounter-sample
r = 0.3;                               % collision radius of each person
NTTC = 30;                              % Number of TTC to sample at first evasive action for each encounter
est_ttc = 1;                             % tells the algorithm whether or not you want to estimate ttc distribution for each encounter
compute_X = 0;
[plots.enc, plots.predpath, plots.pause] =  deal(1, 0, 2^-5); % plot settings

plot_enc = 1;                            % set to one if plots of encounters are wanted
plot_pred_path = 0;                      % enable to simulate each predicted path
disable_crash = 1;                       % disable collision and EA mode to record free movement patterns
pause_length = 2^-5;
use_history = 0;
use_dp_model = 0;                        % set to 1 if you want to use estimated reaction model

for kk = 1:n
kk %#ok<NOPTS>

% parameters for the gamma distribution for the initial steps
speed_mean_0 = 0.16;
speed_var_0  = 0.0015;                     % variance of initial speed
theta_var_0  = 0.015;                      % variance of initial theta

% gamma parametes for initial step as fcn of mean and variance
a_init = speed_mean_0^2/speed_var_0;
b_init = speed_var_0/speed_mean_0;
xinit = 6;                                  % determines how far apart they are at the start
sigma_var0 = 0.3;                           % variance of initial starting position


% variance of the step-size and step-angle
var_step0 = 0.001;
var_theta0 = 0.03;                           % determins the amount of random variability of the step-angle at each iteration.

% EA behaviour modification parameters
thetamod_k = 0.4;              % determines how quickly the person can change direction after detecting the other driver.
thetamod_s = 1.6;              % makes modification of behaviour more dramatic for very small ydiff, but less dramatic for large ydiff
thetamod_p = 3;                % amplifies difference in reaction between small and large ydiff


% variables that collects data about encounters
dist_fea = inf*ones(N,1);               % saves danger index at time of first evasive action
ttc_fea = inf*ones(N,1);                % saves TTC at FEA
stoch_ttc_fea = nan*ones(N, NTTC);      % saves estimates of TTC distributions; each row contains NTTC simulated ttc values
dist_min_ea = inf*ones(N,1000);         % saves min sep. dist. over EA frames
ttc_min_ea = inf*ones(N,1000);          % saves min TTC over EA frames
dist_min_nea = inf*ones(1,N);           % saves min sep. dist. for encounters of type NEA
dist_min = inf*ones(1,N);               % saves min sep. dist. for all encounters
ttc_min = inf*ones(1,N);                % saves min TTC for all encounters
enc_status = ones(0);                   % saves interaction status of each timeframe
min_ttcvec = ones(0);                   
min_distvec = ones(0);
distvec = ones(0);

enc_type = zeros(1,N);                  % vector describing type of each encounter
min_speed = 0.01;                       % minimum stepsize/speed
stability_fac = 0.06;                   % parameter that determines how strongly course gets corrected once disrupted

react_A_save = ones(N,500)*Inf;
react_B_save = ones(N,500)*Inf;
A_save = ones(N,500)*Inf;
B_save = ones(N,500)*Inf;
A_theta_save = ones(N,500)*Inf;
B_theta_save = ones(N,500)*Inf;
A_stepsize_save = ones(N,1);
B_stepsize_save = ones(N,1);

if disable_crash==1
    save_pos_A = inf*ones(N,500);
    save_step_A = inf*ones(N,500);
    save_theta_A = inf*ones(N,500);
    save_dtheta_A = inf*ones(N,500);
    save_dtheta_A(:,1) = 0;
    save_pos_B = inf*ones(N,500);
    save_step_B = inf*ones(N,500);
    save_theta_B = inf*ones(N,500);
end

aa =1;
bb =1;

for i=1:N
decision_state = 0;
i
%%%%%%%%%%%%% initiation of encounter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn = 0;
EA_index = 1;                                                   % tracks number of EA frame
first_detec = 0;                                                % 0 --> neither have not detected yet
speed = min_speed*[1 1] + gamrnd(a_init*[1 1], b_init);   % determines average speed of each vehicle
speed0 = speed;
theta = normrnd([0,pi], theta_var_0);
dtheta = [0,0];
% detection parameters
farsight = 40;      % determines how far drivers extrapolate into future to determine if EA is necessary
shift_d = 2*0.3*2.5;% determines the inflection point of the dist. factor
pow_d = 4;          % determines how quickly dist. factor goes to 1 as min. dist. goes to 0
shift_t = 15;       % determines the inflectopn point of the ttc factor
pow_t = 2;          % determines how quickly ttc factor goes to 1 as ttc goes to 0
detec_amp = 0.45;
detec_par = [shift_d, pow_d, shift_t, pow_t, detec_amp];

tolerance = 0.4*[0.16*0.1, 0.03, 0.08];

%%%% driver properties %%%%
% determine inertia of each driver; how twitchy can they be in their
% movements?
var_theta = var_theta0*gamrnd(1^2/0.04, 0.04/1); 
var_step = var_step0*gamrnd(1^2/0.02, 0.02/1);
d1_theta = +0.5-4.0*gamrnd(1^2/0.1, 0.1/1);           % determines where driver A starts turning
p1_theta = -1.0*gamrnd(1^2/0.1, 0.1/1);               % determines how quickly driver A makes the turn
amp1_theta = pi/2;                                    % determins the angle driver A expects to have after the turn
d1_step = pi/4*0.6*gamrnd(1^2/0.02, 0.02/1);          % angle for which maximum deccelaration occurs
p1_step = 15*gamrnd(1^2/0.2, 0.2/1);              
amp1_step = 0.06*gamrnd(1^2/0.15, 0.15/1);            % determines maximum amount of deceleration during turn

alertness = betarnd(40,3,[1,2]);                      % determines how likely each driver is to notice threat

theta_mod_par = [d1_theta, p1_theta, amp1_theta];      % parameters of s-function that computes exp_theta(A0)
speed_mod_par = [d1_step, p1_step, amp1_step];         % parameters of exp-function that computes exp_step(A0,theta)

% where does A start and end turn
x0 = normrnd(-4.8,0.04);
par = logpar(5,0.5);    
x1 =  x0 + lognrnd(par(1),par(2));
% collect properties characterizing each driver into singe cell array
driver_prop = {speed0, var_step, var_theta, ...
               stability_fac, min_speed, theta_mod_par, speed_mod_par};

% parameters that determine momentum.
weight = gam_rnd(1,0.05,2);
x = struct('a',{10,20});

max_delta(1).dspeed = beta_rnd(5e-4*weight(1), 8e-4);
max_delta(2).dspeed = beta_rnd(5e-4*weight(2), 8e-4);
max_delta(1).dtheta = beta_rnd(0.05*weight(1), 8e-4);
max_delta(2).dtheta = beta_rnd(0.05*weight(2), 8e-4);
max_delta(1).d2theta = deg2rad(4);
max_delta(2).d2theta = deg2rad(4);

%%% simulate initial state
A_real = -xinit;
B_real = xinit;
A_im = normrnd(-2.5,sigma_var0);
B_im = normrnd(0,sigma_var0);
A0 = A_real + 1i*A_im;
B0 = B_real + 1i*B_im;

% collect information about states of drivers
[S(1).pos,S(1).theta,S(1).speed,S(1).dtheta] = deal(A0, theta(1), speed0(1), dtheta(1)); 
[S(2).pos,S(2).theta,S(2).speed,S(2).dtheta] = deal(B0, theta(2), speed0(2), dtheta(2));

%%% data collection variables
% interaction status
int_state.a = 0;   
int_state.b = 0;
encounter_classifier = 1; % tracks status: sign indicates interaction status 
                          % + --> no interaction, - --> interaction); 1 --> no collision, 2 --> collision
danger_enc_i = [];        % saves danger index associated with each timestep
ttc_enc_i = [];           % saves danger index associated with each timestep

% if disable_crash==1
%   engage_EA = double(imag(A0)<4 || real(B0)>-xinit); 
% else
%   engage_EA = double(real(A0) < xinit   && imag(A0)<4 && real(B0) > -xinit && imag(A0)<imag(B0)+1.5 && real(A0)<real(B0)+1.5);
% end
%     while  imag(A0)<4 || real(B0)>-xinit 
    while continue_simulation(A0,B0,xinit,disable_crash)==1
        nn = nn + 1;
        D = norm(A0-B0);
        if disable_crash==1
            save_pos_A(i,nn) = A0;
            save_step_A(i,nn) = speed(1);
            save_theta_A(i,nn) = theta(1);
            save_pos_B(i,nn) = B0;
            save_step_B(i,nn) = speed(2);
            save_theta_B(i,nn) = theta(2);
        end
        
%       dp = p_engage_EA(A0,B0,stepsize,theta, detec_par);
        if int_state.a==0 || int_state.b==0
            
            [dp,t_min,d_min] = pd_separate_reaction(A0,B0,speed0,theta, driver_prop, detec_par, farsight);
            
            if use_dp_model == 1
                predictor = [1,t_min,d_min];
                dp = 1./(1 + exp(-predictor*b));
            end
            
            int_state.a = max(double(rand(1)<dp*alertness(1)),int_state.a);
            int_state.b = max(double(rand(1)<dp*alertness(2)),int_state.b);
            enc_status(end+1) = max(int_state.a,int_state.b);
            min_ttcvec(end+1) = t_min;
            min_distvec(end+1) = d_min;
            detec_stat = max(int_state.a,int_state.b);
            distvec(end+1) = D;
        end
        
        %enc_status(i,nn) = react;
        %dp = 0.98*exp(-0.4*norm(A0-B0));
        if disable_crash == 1 
            dp=0; 
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if          "drivers are in EA mode"    OR     "EA is engaged by dp roll"
        % and         "they have not yet passed each other"
        if  detec_stat==1 && real(A0) < real(B0) && imag(A0)<imag(B0)+0.5 && disable_crash==0
            % saving information in vectors, and updating status
            nn
            if first_detec == 0 && est_ttc == 1 % a loop that collects information at first moment of evasive action
                
                if use_history==1
                    [row, col] = find_sim_sit(A0, speed(1),theta(1), historyA, tolerance);
                else
                    row=0;
                end
                
                for k=1:NTTC 
                    %pred_paths = inf*ones(NTTC,500);
                    if k>length(row) || use_history==0 %number of observed walks satisfying init. cond.
                        ttc = ttc_simulator_double_momentumV2(A0,B0, speed, theta, r, driver_prop, plot_pred_path);
                        stoch_ttc_fea(i,k) = ttc;
%                         ttc
%                         pause(.1)
                    else
                        ttc = ttc_simulator_history(A0,B0,speed,theta,r,driver_prop,historyA{1},row,col,k,plot_pred_path);
                        stoch_ttc_fea(i,k) = ttc;
%                         ttc;
%                         pause(1)
                        %                     ttc = ttc_simulator_double_momentumV2(A0,B0, stepsize, theta, r, driver_prop, plot_setting);
                    end
                end
                ttc_fea(i) = calc_ttc(S,r);
                dist_fea(i) = norm(A0-B0) - 2*r;
                dist_min_EA(i,EA_index) = norm(A0-B0) - 2*r;
                ttc_min_ea(i,EA_index) = ttc_fea(i);
                EA_index = EA_index + 1;
                
            end
            
            first_detec = 1; % set to 1 so that above loop only runs on first moment of detection
            int_state.a = 1;
            encounter_classifier = -1;
            danger_index = D - 2*r;
            danger_enc_i(length(danger_enc_i) + 1) = danger_index;
            ttc_enc_i(length(ttc_enc_i) + 1) = calc_ttc(S,r);
            
            % take next step
%            if detec_stat_A==1 && detec_stat_B==1
                newstep =  evasive_step_coulombs_law(A0,B0,speed,theta, 5, 5);
                A1 = newstep{1}(1);
                B1 = newstep{1}(2);
                theta = newstep{3};
               

%            elseif detec_stat_A==1 && detec_stat_B==0
%                 newstep =  evasive_step_coulombs_law(A0,A_save(i,nn-1),B0,B_save(i,nn-1),stepsize, 1, 1);
%                 decision_state = newstep{2};
%                 A1 = newstep{1}(1);    
%                 B1 = newstep{1}(2);
%                 stepsize = newstep{2};
%                 theta = newstep{3};
%                 
%             elseif detec_stat_A==0 && detec_stat_B==1
 
                
            dist_min_EA(i,EA_index) = norm(A1-B1) - 2*r;
            ttc_min_ea(i,EA_index) = compute_ttc(A1,B1,speed,theta,r);
            EA_index = EA_index + 1;

            plot_pos(S, plots, xinit, r, detec_stat_A, detec_stat_B) 
            
            if norm(A1-B1) < 2*r && disable_crash==0 % if satisfied, then colliosion with evasive action has ocurred
                hold on
                title(sprintf(" decision_state = %d", decision_state))
                hold off
                danger_index = norm(A1-B1) - 2*r;
                dist_min(i) = danger_index;
                ttc_min(i) = calc_ttc(S,r);
                danger_enc_i(length(danger_enc_i) + 1) = danger_index;
                ttc_enc_i(length(ttc_enc_i) + 1) = ttc_min(i);
                encounter_classifier = -2;             % set encounter to type crash with attempted evasive action
                break
            end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% no detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            % saving information in vectors, and updating status

            danger_index = D - 2*r;
            danger_enc_i(length(danger_enc_i) + 1) = danger_index;
            ttc_enc_i(length(ttc_enc_i) + 1) = calc_ttc(S,r);
            
            % find desired angle and speed
            E_thetaA = normrnd(expected_angle_A(S, x0, x1), var_theta0);
            % add small value to argument to make cars slow down before
            % they start to make turn
            E_speedA = normrnd(expected_speed_A(real(A0) + .3, x0, x1, speed0(1)), var_step0);
            E_thetaB = normrnd(pi,0.1);
            E_speedB = normrnd(speed0(2),0.01);
            
            % take step and update states
            S(1) = take_step(S(1), E_speedA, E_thetaA, max_delta(1)); 
            S(2) = take_step(S(2), E_speedB, E_thetaB, max_delta(2));
            
            A1          = S(1).pos;
            speed(1)    = S(1).speed;
            dtheta(1)   = S(1).dtheta;
            theta(1)    = S(1).theta;
            
            B1          = S(2).pos;
            speed(2)    = S(2).speed;
            dtheta(2)   = S(2).dtheta;
            theta(2)    = S(2).theta;
            
            if min(S(1).speed,S(2).speed)<=1e-4
                break
            end
            
            D = norm(S(1).pos - S(2).pos);
            ttc = calc_ttc(S,r);
            
            plot_pos(S, plots, xinit, r, int_state)
            
            if D < 2*r && disable_crash==0% if satisfied, collision has occured before evasive action is taken
                hold on; title("collision"); hold off
                danger_index = D - 2*r;

                %%% compute ttc (which will now be negative) and danger_FEA
                ttc = ttc_simulator_double_momentum_improved(A1,B1,speed,theta,min_speed,var_step,r, var_theta, stability_fac);

                dist_fea(i) = danger_index;
                ttc_fea(i) = ttc;
                stoch_ttc_fea(i,:) = ttc*ones(1,NTTC);
                danger_enc_i(length(danger_enc_i) + 1) = danger_index; %#ok<*SAGROW>
                ttc_enc_i(length(ttc_enc_i) + 1) = compute_ttc(A1,B1,speed,theta,r);
                dist_min(i) = danger_index;
                ttc_min(i) = compute_ttc(A1,B1,speed,theta,r);
                encounter_classifier = 2;                                  % indicates a collision with no evasive action attempted
                break
            end
        end
        react_A_save(i,nn)  = int_state.a;
        react_B_save(i,nn)  = int_state.b;
        A_save(i,nn) = A1;
        B_save(i,nn) = B1;
        A_theta_save(i,nn)  = theta(1);
        B_theta_save(i,nn)  = theta(2);
        if nn>1
            A_dtheta_save(i,nn) = A_theta_save(i,nn) - A_theta_save(i,nn-1);
        end

        A0 = A1;
        B0 = B1;
    end
    
    A_stepsize_save(i) = speed(1);
    B_stepsize_save(i) = speed(2);

    enc_type(i) = encounter_classifier;
    dist_min(i) = min(danger_enc_i);
    ttc_min(i) = min(ttc_enc_i);

    if encounter_classifier == 1  % no evasive action and no collision has occured
        dist_min_nea(i) = min(danger_enc_i);
    end
    sum(enc_type==2);

end

if disable_crash == 1
    historyA{1} = save_pos_A;
    historyA{2} = save_step_A;
    historyA{3} = save_theta_A;
end


%plot_encounter(A_save, B_save,react_A_save, react_B_save, find(enc_type==-2), .15, xinit, r)




sum(A_save(enc_id(i),:)<inf)
























sum(enc_type==-2)
sum(enc_type==2)
% remove all elements/rows corresponding to encounters with no EA
dist_fea(enc_type==1,:) =[];
X = stoch_ttc_fea;
stoch_ttc_fea(enc_type==1,:) = [];
dist_min_EA(enc_type>0,:) = [];
dist_min_EA = min(dist_min_EA,[],2);
ttc_min_ea( enc_type > 0 ,:) = [];
ttc_min_ea = min(ttc_min_ea,[],2);
dist_min_nea = dist_min_nea( enc_type > 0 );

%%%%%%% Compute minimum ttc for all encounters where there was no evasive action
if est_ttc == 1 && compute_X ==1

N_enc_type1 = sum(enc_type==1);
A_NEA = A_save(enc_type==1,:);
B_NEA = B_save(enc_type==1,:);
diff_NEA = A_NEA - B_NEA;
real_diff_NEA = real(diff_NEA);
norm_diff_NEA = sqrt(conj(diff_NEA).*diff_NEA);
ttc_dist_save = ones(N_enc_type1, NTTC)*Inf;

for i=1:N_enc_type1
    %i
    ttc_dist = ones(500, NTTC)*Inf;
    ttc_minimum = Inf;
    min_ttc_index = 1;
%     if min(norm_diff_NEA(i,:))>2
%         continue
%     end

    for j=1:100
        t = norm_diff_NEA(i,j)/(A_stepsize_save(i) + B_stepsize_save(i));
        if real_diff_NEA(i,j)<0   &   t < 6    &     min(norm_diff_NEA(i,:)) < 4
                A0 = A_NEA(i,j);
                B0 = B_NEA(i,j);
                theta = [A_theta_save(i,j), B_theta_save(i,j)];
                speed = [A_stepsize_save(i), B_stepsize_save(i)];

                for k=1:NTTC
                    ttc = ttc_simulator_double_momentum_improved(A0,B0,speed,theta,min_speed,var_step*aa,r, var_theta*bb, stability_fac,plot_sim_walks);
                    if ttc< Inf
                        pause(0.0)
                    end
                    if ttc < ttc_minimum
                        ttc_minum=ttc;
                        min_ttc_index = j;    % keeps track of which moment is most dangerous
                    end
                    ttc_dist(j,k) = ttc;
                    if ttc<0
                        TTC = 0;
                        AA = A0;
                        BB = B0;
                        step0 = speed;
                        theta0 = theta;
                    end

                end



    end
    ttc_dist_save(i,:) = ttc_dist(min_ttc_index,:); % find ttc-distribution from moment with lowest minimum ttc.
    end
end

X(enc_type==1,:) = ttc_dist_save;   % vector that contains ttc-data from all encounters.
end

avg_cond_ttc = mean(stoch_ttc_fea,2,'omitnan');
p_hypo_coll = sum(1-isnan(stoch_ttc_fea'))'/NTTC;
weighted_ttc_fea = avg_cond_ttc.*(1+exp(-3*p_hypo_coll));
% ttc and stochastic-ttc measurements
all_data{kk,1} = stoch_ttc_fea;
all_data{kk,2} = weighted_ttc_fea;
all_data{kk,3} = ttc_fea;
all_data{kk,4} = ttc_min';
all_data{kk,5} = ttc_min_ea;

% minim distance measurements
all_data{kk,6} = dist_fea;
all_data{kk,7} = dist_min_EA;
all_data{kk,8} = dist_min';

% other information
all_data{kk,9} = sum(enc_type==2);                                               % saves p_nea_coll
all_data{kk,10} = sum(enc_type==-2);                                             % saves p_ea_coll
all_data{kk,11} = (sum(enc_type==-1) + sum(enc_type==-2) + sum(enc_type==2))/N;  % saves p_interactive
all_data{kk,12} = (sum(enc_type==-1)+sum(enc_type==-2))/N;                       % saves p_ea


end