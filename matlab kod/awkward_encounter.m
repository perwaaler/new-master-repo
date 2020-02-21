
clear all
pause_length = 1/4/2/2/2/2/2/2;
plotting = 1;
% parameters for distribution of step-angle
m = 5; % increase m --> less variance in stepsize
alpha0 = m;
beta0 = m;

% parameters for distribution of stepsize. EX = a*b, VX = a*b^2
a = 0.3*2;
b = 0.2/2;

% parameters relating to detection probability, pd = exp(-s*(D + d)^p)
d = .60;
s = 0.35;   % decrease to make it harder to detect
p = 1.3;   % increase to make it unlikelier to detect for large d, but more likely for small d
plot(linspace(0,10,100),exp(-s*(linspace(0,10,100)+d).^p))

% Collision Avoidance Parameters: Parameter-Modification Factor = exp( -s*((x + d)/k)^-p )
d_ca1 = 0.1;        % increase to reduce PMF
k_ca1 = 1.4;    % increase to decrease modification for x<k-d and reduce modification for x>k-d
p_ca1 = 1.6;    % increases effect of avoidance modifier when 
s_ca1 = 0.2

d_ca2 = 0.1;        % increase to reduce PMF
k_ca2 = 2.6;    % increase to decrease modification for x<k-d and reduce modification for x>k-d
p_ca2 = 1.8;    % increases effect of avoidance modifier when 
s_ca2 = 2

theta_range = pi/4; % range of the distribution of stepangles

r = 0.5; % collision radius of each person
xinit = 10;
sigma = 0.3;
N=50;
danger_FEA = zeros(1,N);         % danger index at time of first evasive action
danger_CR = zeros(1,N);          % danger index at time of conflict resolution
danger_nodetec = zeros(1,N); % highest index of danger for non-detection encounters
enc_type = zeros(1,N);      % vector containing encounter type of each encounter
for i=1:N
% initiation
A_real = -xinit;
B_real = xinit;
A_im = normrnd(0,sigma);
B_im = normrnd(0,sigma);
A = A_real+1i*A_im;
B = B_real+1i*B_im;
% data collection variables
detection_status = 0; 
encounter_classifier = 1; % 1 --> (no detection, no crash), -1 --> (detection, no crash) 2 --> (no detection, crash), -2 --> (detection, crash)
danger_vec = []; % collects the danger index associated with each timestep
detec_vec = [];    % interaction vector, 1 --> detection at timeframe, 0 --> no detection at timeframe

    while real(A) < xinit && real(B) > -xinit
        D = norm(A-B);
        pd = exp(-s*(D+d)^p);
        %%%%%%%%%%%%%%%% detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  (detection_status==1 | rand(1) < pd) && real(A) < real(B)
            % saving information in vectors, and updating status
            detection_status = 1;
            encounter_classifier = -1;
            detec_vec(length(detec_vec) + 1) = 1;
            danger_index = D - 2*r;
            danger_vec(length(danger_vec) + 1) = danger_index;

            % take next step
            pred_posA = A + a*b;
            pred_posB = B - a*b;
            pred_diff = pred_posA - pred_posB;
            pred_dist = norm(pred_diff);
            alpha_modifier = mod_fac(pred_dist,d_ca1,k_ca1,p_ca1,s_ca1);
            theta_modifier = mod_fac(abs(imag(pred_diff)),d_ca2,k_ca2,p_ca2,s_ca2);
            stepsize = gamrnd(a*alpha_modifier*[1 1],b*[1 1]); % step-size
            theta = betarnd(alpha0*theta_modifier*[1 1], beta0*[1 1]);
            theta = sign(imag(pred_diff))*[-1 -1].*(theta*theta_range - 0.5*theta_range); % step-angle
            A = A + stepsize(1)*exp(1i*theta(1));
            B = B - stepsize(2)*exp(1i*theta(2));

            % plotting %
            if plotting == 1
                xlim([-xinit,xinit])
                ylim([-4,4])
                plot([r*cos(linspace(0,2*pi,50))+real(A)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A)])
                title("detection")
                hold on
                plot([r*cos(linspace(0,2*pi,50))+real(B)]+1i*[r*sin(linspace(0,2*pi,50))+imag(B)])
                xlim([-xinit,xinit])
                ylim([-4,4])
                hold off
                pause(pause_length)
            end
            if norm(A-B) < 2*r
                hold on
                title("collision")
                hold off
                detec_vec(length(detec_vec) + 1) = 1;
                danger_index = norm(A-B) - 2*r;
                danger_vec(length(danger_vec) + 1) = danger_index;
                encounter_classifier = -2;
                break
            end
        %%%%%%%%%%%%%%%% no detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            % saving information in vetctors, and updating status
            detec_vec(length(detec_vec) + 1) = 0;
            danger_index = D - 2*r;
            danger_vec(length(danger_vec) + 1) = danger_index;

            % taking next step
            theta = betarnd(alpha0*[1 1],beta0*[1 1]);
            theta = theta*theta_range - 0.5*theta_range; % step-angle
            stepsize = gamrnd(a*[1 1],b*[1 1]); % step-size
            A = A + stepsize(1)*exp(1i*theta(1));
            B = B - stepsize(2)*exp(1i*theta(2));
            D = norm(A-B);

            % plotting %
            if plotting == 1
                xlim([-xinit,xinit])
                ylim([-4,4])
                plot([r*cos(linspace(0,2*pi,50))+real(A)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A)])
                title("no detection")
                hold on
                plot([r*cos(linspace(0,2*pi,50))+real(B)]+1i*[r*sin(linspace(0,2*pi,50))+imag(B)])
                xlim([-xinit,xinit])
                ylim([-4,4])
                hold off
                pause(pause_length)
            end
            if norm(A-B) < 2*r 
                hold on
                title("collision")
                hold off
                detec_vec(length(detec_vec) + 1) = 1;
                danger_index = norm(A-B) - 2*r;
                danger_vec(length(danger_vec) + 1) = danger_index; 
                encounter_classifier = 2;
                break
            end
        end
    end

    enc_type(i) = encounter_classifier;
    detec_dangerind = danger_vec(find(detec_vec)); % contains danger index for all times with detection status
    if encounter_classifier == 1
        danger_nodetec(i) = min(danger_vec);
    else
        danger_FEA(i) = detec_dangerind(1);
        danger_CR(i) = detec_dangerind(end);    
    end

end
%%
save('danger_CR2million.mat','danger_CR')
save('danger_FEA2million.mat','danger_FEA')
%% parameters used when doing 10^6 iterations.

% parameters for distribution of step-angle
m = 5; % increase m --> less variance in stepsize
alpha0 = m;
beta0 = m;

% parameters for distribution of stepsize. EX = a*b, VX = a*b^2
a = 0.3*2;
b = 0.2/2;

% parameters relating to detection probability, pd = exp(-s*(D + d)^p)
d = .60;
s = 0.35;   % decrease to make it harder to detect
p = 1.3;   % increase to make it unlikelier to detect for large d, but more likely for small d
plot(linspace(0,10,100),exp(-s*(linspace(0,10,100)+d).^p))

% Collision Avoidance Parameters: Parameter-Modification Factor = exp( -s*((x + d)/k)^-p )
d_ca1 = 0.1;        % increase to reduce PMF
k_ca1 = 1.4;    % increase to decrease modification for x<k-d and reduce modification for x>k-d
p_ca1 = 1.6;    % increases effect of avoidance modifier when 
s_ca1 = 0.2

d_ca2 = 0.1;        % increase to reduce PMF
k_ca2 = 2.6;    % increase to decrease modification for x<k-d and reduce modification for x>k-d
p_ca2 = 1.8;    % increases effect of avoidance modifier when 
s_ca2 = 2

theta_range = pi/4; % range of the distribution of stepangles

r = 0.2; % collision radius of each person
xinit = 10;
sigma = 0.3;
N=1000000;
danger_FEA = zeros(1,N);         % danger index at time of first evasive action
danger_CR = zeros(1,N);          % danger index at time of conflict resolution
danger_nodetec = zeros(1,N); % highest index of danger for non-detection encounters
enc_type = zeros(1,N);      % vector containing encounter type of each encounter
%% Data from experiment with 10^6 encounters along with confidence intervals and estimates.
N = 10^6;
n_coll = 110;
p_coll = 110/10^6;
sdev = sqrt( p_coll*(1-p_coll)/N )
ci = p_coll + [-1,1]*1.96*sdev

pFEA = sum(danger_FEA<0)/length(danger_FEA)
pCR = sum(danger_CR<0)/length(danger_CR)
%% Data from experiment with 2*10^6 encounters along with confidence intervals
N = 2*10^6;
n_coll = 249;
p_coll = n_coll/N;
sdev = sqrt( p_coll*(1 - p_coll)/N )
ci = p_coll + [-1,1]*1.96*sdev




