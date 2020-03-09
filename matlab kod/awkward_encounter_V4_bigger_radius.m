% simulation of encounters between two vehicles.

n = 1;                                % number of encounter-samples desired
all_data = cell(n, 12);                 % collects data from each encounter-sample
N = 10000;                                % number of encounters
r = 0.3;                                % collision radius of each person
NTTC = 100;                             % Number of TTC to sample at first evasive action for each encounter
est_ttc = 0;                            % tells the algorithm whether or not you want to estimate ttc distribution for each encounter
compute_X = 0;
plotting = 0;                           % set to one if plots of encounters are wanted
sausage = 0;
plot_sim_walks = 0;

for kk = 1:n
kk %#ok<NOPTS>

pause_length = 2^-5;

% parameters for distribution of stepsize. EX = a*b, VX = a*b^2
step_par = [0.30*2 0.2/2];

% parameters for the gamma distribution for the initial steps
step_mean_init = 0.125*0.93;
step_var_init = 0.003;
theta_var_init = 0.015;                      % variance of initial theta


% parametes as fcn of mean and variance
a_init = step_mean_init^2/step_var_init;
b_init = step_var_init/step_mean_init;

% variance of the step-size
var_step = 0.001*2;
var_theta = 0.03*2;                           % determins the amount of random variability of the step-angle at each iteration.

% parameters relating to detection probability dp; dp = A*exp(-s*(D + d)^p)
% where D is distance between vehicles
Amp = 0.93;      % amplitude; reduce to make dp smaller everywhere
d = 0;           % determines where dp is made smaller and larger by p; dp gets larger for D smaller than d, and smaller for D greater than d
s = 0.26;        % increase to make detection less likely
p = 2;           % makes dp larger r smaller depending on whether D>d or D<d. See above comment.
% plot(linspace(0,10,100),Amp*exp(-(s * linspace(0,10,100) ).^p))

thetamod_k = 0.3;                       % determines how quickly the person can change direction after detecting the other person
thetamod_s = 1.6;                       % makes modification of behaviour more dramatic for very small ydiff, but less dramatic for large ydiff
thetamod_p = 3;                         % amplifies difference in reaction between small and large ydiff
xinit = 6;                              % determines how far apart they are at the start
sigma_var0 = 0.3;                            % variance of initial starting position

dist_FEA = inf*ones(N,1);            % danger index at time of first evasive action
ttc_FEA = inf*ones(N,1);
stoch_ttc_FEA = inf*ones(N, NTTC);              % Saves danger at first moment of evasive action; each row contains NTTC sampled ttc values
dist_min_EA = inf*ones(N,1000);       % collects danger measurements at each time frame where there is evasive action
ttc_min_EA = inf*ones(N,1000);       % collects danger measurements at each time frame where there is evasive action
dist_min_NEA = inf*ones(1,N);    % highest index of danger for non-detection encounters (i.e. type 1 or 2 encounters)
dist_min = inf*ones(1,N);              % collects data from most dangerous timeframe from each encounter
ttc_min = inf*ones(1,N);

enc_type = zeros(1,N);                  % vector containing encounter type of each encounter
min_stepsize = 0.02;                    % set to make it impossible for a stepsize to be smaller than this
stability_fac = 0.04;

A_save = ones(N,500)*Inf;
B_save = ones(N,500)*Inf;
A_theta_save = ones(N,500)*Inf;
B_theta_save = ones(N,500)*Inf;
A_stepsize_save = ones(N,1);
B_stepsize_save = ones(N,1);
aa =1;
bb =1;

for i=1:N
i
%%% initiation of encounter
nn = 0;
EA_index = 1;                                                         % variable used to find most dangerous moment during attempt to avoid collision
first_detection = 0;                                                  % They have not yet detected eachother
stepsize = min_stepsize*[1 1] + [gamrnd(a_init*[1 1],b_init)];        % determines average speed of each vehicle
theta = normrnd([0,0], theta_var_init);

%%% initial positions
A_real = -xinit;
B_real = xinit;
A_im = normrnd(0,sigma_var0);
B_im = normrnd(0,sigma_var0);
A0 = A_real + 1i*A_im;
B0 = B_real + 1i*B_im;

%%% data collection variables
detection_status = 0;
encounter_classifier = 1; % tracks status: sign indicates interaction status (+ --> no interaction, - --> interaction), number indicates collision status (1 --> no collision, 2 --> collision)
danger_enc_i = [];        % saves danger index associated with each timestep
ttc_enc_i = [];        % saves danger index associated with each timestep
counter = 0;
    while real(A0) < xinit && real(B0) > -xinit
        counter = counter + 1;
        nn = nn + 1;

        D = norm(A0-B0);
        Dim = norm(imag(A0-B0));                  % vertical part of the distance; if large, evasive action should be unlikely
        dp = Amp*exp(-(s*D)^p)*exp(-(Dim/0.8000)^3); % detection probability

%%%%%%%%%%%%%%%%%%%%%%%% detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  (detection_status==1 | rand(1) < dp) && real(A0) < real(B0)
            % saving information in vectors, and updating status
            if first_detection == 0 & est_ttc == 1 % a loop that collects information at first moment of evasive action
                for k=1:NTTC
                    ttc = ttc_simulator_double_momentum_improved(A0,B0,stepsize,theta,min_stepsize,var_step*aa,r, var_theta*bb, stability_fac,plot_sim_walks);
                    stoch_ttc_FEA(i,k) = ttc;
                end
                ttc_FEA(i) = compute_ttc(A0,B0,stepsize,theta,r);
                if ttc_FEA(i) < min(stoch_ttc_FEA(i,:))
                    sausage = sausage+1
                    ttc_FEA(i)
                    min(stoch_ttc_FEA(i,:))
                    pause(1)
                                    xlim([-xinit,xinit])
                    ylim([-4,4])
                    Apast = A0 - stepsize(1)*exp(1i*theta(1));
                    Bpast = B0 + stepsize(2)*exp(1i*theta(2));
                    Afut = A0 + stepsize(1)*exp(1i*theta(1))*ttc_FEA(i);
                    Bfut = B0 - stepsize(2)*exp(1i*theta(2))*ttc_FEA(i);
                    clf
                    plot([r*cos(linspace(0,2*pi,50))+real(Apast)]+1i*[r*sin(linspace(0,2*pi,50))+imag(Apast)],'b');
                    hold on
                    plot([r*cos(linspace(0,2*pi,50))+real(A0)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A0)],'b');
                    plot([r*cos(linspace(0,2*pi,50))+real(Afut)]+1i*[r*sin(linspace(0,2*pi,50))+imag(Afut)],'r');

                    title("detection")
                    hold on
                    plot([r*cos(linspace(0,2*pi,50))+real(B0)]+1i*[r*sin(linspace(0,2*pi,50))+imag(B0)],'g')
                    plot([r*cos(linspace(0,2*pi,50))+real(Bpast)]+1i*[r*sin(linspace(0,2*pi,50))+imag(Bpast)],'g')
                    plot([r*cos(linspace(0,2*pi,50))+real(Bfut)]+1i*[r*sin(linspace(0,2*pi,50))+imag(Bfut)],'black')

                    xlim([-xinit,xinit])
                    ylim([-4,4])
                    hold off
                    %saveas(gcd,sprintf('fig_%d_%d.jpg',i,counter))
                    %savefig(sprintf('fig_%d_%d',i,counter))
                    pause(4)

                    %ttc_simulator_double_momentum_improved(A0,B0,stepsize,theta,min_stepsize,var_step*aa,r, thetavar*bb, stability_fac,1)
                end

                pause(0)

                dist_FEA(i) = norm(A0-B0) - 2*r;
                dist_min_EA(i,EA_index) = norm(A0-B0) - 2*r;
                ttc_min_EA(i,EA_index) = compute_ttc(A0,B0,stepsize,theta,r);
                EA_index = EA_index + 1;
            end

            first_detection = 1;                            % set to 1 so that above loop only runs on first moment of detection
            detection_status = 1;
            encounter_classifier = -1;
            danger_index = D - 2*r;
            danger_enc_i(length(danger_enc_i) + 1) = danger_index;
            ttc_enc_i(length(ttc_enc_i) + 1) = compute_ttc(A0,B0,stepsize,theta,r);

            % take next step
            [A1, B1, theta] = take_evasive_step(A0,B0,step_par,stepsize,var_step, min_stepsize,theta,var_theta,stability_fac,thetamod_k,thetamod_p,thetamod_s);

            dist_min_EA(i,EA_index) = norm(A1-B1) - 2*r;
            ttc_min_EA(i,EA_index) = compute_ttc(A1,B1,stepsize,theta,r);
            EA_index = EA_index + 1;

            if plotting == 1
                xlim([-xinit,xinit])
                ylim([-4,4])
                plot([r*cos(linspace(0,2*pi,50))+real(A1)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A1)])
                title("detection")
                hold on
                plot([r*cos(linspace(0,2*pi,50))+real(B1)]+1i*[r*sin(linspace(0,2*pi,50))+imag(B1)])
                xlim([-xinit,xinit])
                ylim([-4,4])
                hold off
                pause(pause_length)
            end
            if norm(A1-B1) < 2*r                % this part is entered, then a colliosion with evasive action has occured
                hold on
                title("collision")
                hold off
                danger_index = norm(A1-B1) - 2*r;
                dist_min(i) = danger_index;
                ttc_min(i) = compute_ttc(A1,B1,stepsize,theta,r);
                danger_enc_i(length(danger_enc_i) + 1) = danger_index;
                ttc_enc_i(length(ttc_enc_i) + 1) = compute_ttc(A1,B1,stepsize,theta,r);
                encounter_classifier = -2;             % set encounter to type crash with attempted evasive action
                break
            end
%%%%%%%%%%%%%%%%%%%%%%%% no detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            % saving information in vectors, and updating status

            danger_index = D - 2*r;
            danger_enc_i(length(danger_enc_i) + 1) = danger_index;
            ttc_enc_i(length(ttc_enc_i) + 1) = compute_ttc(A0,B0,stepsize,theta,r);

            % taking next step
            theta = theta + normrnd(-stability_fac*theta,var_theta);
            stepsize = stepsize + normrnd([0 0],var_step);
            stepsize(1) = max(min_stepsize,stepsize(1));
            stepsize(2) = max(min_stepsize,stepsize(2));
            A1 = A0 + stepsize(1)*exp(1i*theta(1));
            B1 = B0 - stepsize(2)*exp(1i*theta(2));
            D = norm(A1-B1);
            ttc = compute_ttc(A1,B1,stepsize,theta,r);
            if plotting == 1
                xlim([-xinit,xinit])
                ylim([-4,4])
                plot([r*cos(linspace(0,2*pi,50))+real(A1)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A1)])
                title("no detection")
                hold on
                plot([r*cos(linspace(0,2*pi,50))+real(B1)]+1i*[r*sin(linspace(0,2*pi,50))+imag(B1)])
                xlim([-xinit,xinit])
                ylim([-4,4])
                hold off
                %saveas(gcf,sprintf('fig_%d_%d.jpg',i,counter))
                %savefig(sprintf('fig_%d_%d',i,counter))
                pause(pause_length)
            end
            if D < 2*r                                           % collision has occured before evasive action is taken
                hold on; title("collision"); hold off
                danger_index = D - 2*r;

                %%% compute ttc (which will now be negative) and danger_FEA
                ttc = ttc_simulator_double_momentum_improved(A1,B1,stepsize,theta,min_stepsize,var_step,r, var_theta, stability_fac);

                dist_FEA(i) = danger_index;
                ttc_FEA(i) = ttc;
                stoch_ttc_FEA(i,:) = ttc*ones(1,NTTC);
                danger_enc_i(length(danger_enc_i) + 1) = danger_index; %#ok<*SAGROW>
                ttc_enc_i(length(ttc_enc_i) + 1) = compute_ttc(A1,B1,stepsize,theta,r);
                dist_min(i) = danger_index;
                ttc_min(i) = compute_ttc(A1,B1,stepsize,theta,r);
                encounter_classifier = 2;                                  % indicates a collision with no evasive action attempted
                break
            end
        end
        A_save(i,nn) = A1;
        B_save(i,nn) = B1;
        A_theta_save(i,nn) = theta(1);
        B_theta_save(i,nn) = theta(2);

        A0 = A1;
        B0 = B1;
    end
    A_stepsize_save(i) = stepsize(1);
    B_stepsize_save(i) = stepsize(2);

    enc_type(i) = encounter_classifier;
    dist_min(i) = min(danger_enc_i);
    ttc_min(i) = min(ttc_enc_i);

    if encounter_classifier == 1                   % i.e. no evasive action and no collision has occured
        dist_min_NEA(i) = min(danger_enc_i);
    end


end
sum(enc_type==-2)

% remove all elements/rows corresponding to encounters with no EA
dist_FEA(enc_type==1,:) =[];
X = stoch_ttc_FEA;
stoch_ttc_FEA(enc_type==1,:) = [];
dist_min_EA(enc_type>0,:) = [];
dist_min_EA = min(dist_min_EA,[],2);
ttc_min_EA( enc_type > 0 ,:) = [];
ttc_min_EA = min(ttc_min_EA,[],2);
dist_min_NEA = dist_min_NEA( enc_type > 0 );

%%%%%%% Compute minimum ttc for all encounters where there was no evasive action
if est_ttc == 1 & compute_X ==1

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
                stepsize = [A_stepsize_save(i), B_stepsize_save(i)];

                for k=1:NTTC
                    ttc = ttc_simulator_double_momentum_improved(A0,B0,stepsize,theta,min_stepsize,var_step*aa,r, var_theta*bb, stability_fac,plot_sim_walks);
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
                        step0 = stepsize;
                        theta0 = theta;
                    end

                end



    end
    ttc_dist_save(i,:) = ttc_dist(min_ttc_index,:); % find ttc-distribution from moment with lowest minimum ttc.
    end
end

X(enc_type==1,:) = ttc_dist_save;   % vector that contains ttc-data from all encounters.
end

% ttc and stochastic-ttc measurements
all_data{kk,1} = stoch_ttc_FEA;
all_data{kk,2} = X;
all_data{kk,3} = ttc_FEA;
all_data{kk,4} = ttc_min';
all_data{kk,5} = ttc_min_EA;

% minim distance measurements
all_data{kk,6} = dist_FEA;
all_data{kk,7} = dist_min_EA;
all_data{kk,8} = dist_min';

% other information
all_data{kk,9} = sum(enc_type==2);                                               % saves p_nea_coll
all_data{kk,10} = sum(enc_type==-2);                                             % saves p_ea_coll
all_data{kk,11} = (sum(enc_type==-1) + sum(enc_type==-2) + sum(enc_type==2))/N;  % saves p_interactive
all_data{kk,12} = (sum(enc_type==-1)+sum(enc_type==-2))/N;                       % saves p_ea

save('data_500enc_r_0_3_extra_rand_safety_level_1','all_data')
 end
