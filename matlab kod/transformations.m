%% Estimation of P(collision)
% this code is used to take simulated data and perform POT method to
% estimate P(collision). Data_matrix  can be either danger_FEA, DAFEA, X,
% or danger_max_EA. data_type=1 --> danger_FEA or DAFEA, data_type=2 --> X,
% data_type=3 --> danger_max_EA
%%
load data_500enc_r_0_3_safety_level_2.mat
%% Initial plots
data_type = 1;
data_matrix = all_data{1,data_type}; % select data. row i should correspond to encounter i, and column j to j'th simulated ttc value (in case of stochastic ttc)

% find encounters with finite ttc values
min_data = data_matrix(min(data_matrix,[],2)<Inf,:);
min_data = min(min_data,[],2);

% find minimum and maximum thresholds based on amount of data to be used

u_minmax = find_threshold(min_data, 0.06, 0.8);

u_l = u_minmax(1);
u_u = u_minmax(2);

clf;
plot(min_data,'.'); hold on 
plot(ones(1,length(min_data))*u_l) 
plot(ones(1,length(min_data))*u_u)
%ylim([-1,100])
title('untransformed data')
%% transforming data and plotting transformed data and thresholds
% transforms
p = 3;

transinv = @(x)1./(3 + x).^p;

transex = @(x)exp(-0.3*(x - 2));

transneg = @(x) -x;

% choose transform to use
trans = @(x) transinv(x);
Nenc = length(data_matrix(:,1));

% transform thresholds
u_min_max = sort(trans([u_l,u_u]));
u_l_trans = u_min_max(1);
u_u_trans = u_min_max(2);
m=5;
U = sort(trans(linspace(u_l, u_u,m)));

trans_data = trans(data_matrix);

clf; plot(trans(min_data),'.'); hold on
plot(1:length(min_data),U'*ones(1,length(min_data)), 'k');plot(ones(1,length(min_data))*trans(0));ylim([trans(30),trans(0)*1.1])
title('transformed data')
legend('transformed thresholds')
%% ensure good initial value!
trans_data = trans_data(:);
exceed = trans_data(trans_data > u_l_trans);
negL = @(par,exceed_data,u) -sum( log(gppdf(exceed_data,par(2),par(1),u)) );
init = fminsearch(@(par) negL(par, exceed, u_l_trans), [3 0.3])
%% Fit models for range of thresholds and show diagnostic plots
shake_guess = 0.1;                 % variance of noise that gets added to initial guess when stuck
Nenc = length(data_matrix(:,1));   % number of encounters
Nexp = length(data_matrix(1,:));   % number of values per row. If surrogate-measure is deterministic, then Nexp=1.
Nbs = 20;                         % number of bootstrapped samples to compute standard error
trans_data = trans(data_matrix);                                                  % number of thresholds used for estimation
                   % vector containing thresholds'

% collects 95% ci's for xi
ci_sigma_u = zeros(2,m)*nan;
ci_xi_u = zeros(2,m)*nan;                                      % collects 95% ci's for xi
ci_p_nea_u = zeros(2,m)*nan;

param_save = zeros(2,m);
pc = zeros(1,m)*nan;                                    % collects estimated collision probability for each threshold
ue_save = zeros(1,m)*nan;                                  % collects estimated upper endpoint
max_data = max(max(trans_data));                           % largest observed value
negL = @(par, exceed_data,u) -sum( log(gppdf(exceed_data,par(2),par(1),u)) );  %negative log likelihood fcn.
logit = 1;                                               % set to plot logarithm of p_nea when magnitude of p_nea varies alot
compute_ci = 1;                    % set equal to one if confidence intervals for xi are desired
qqplot = 1;
pause_length = 0;
% limits for plots of empirical and model distribution functions
xplot_lower = 0;
xplot_upper = 1;
n_eval_cdf = 1000;    % number of points where cdf gets evaluated

for k=1:m
    data = trans_data(:);
    exceed = data(data>U(k));
    param = fminsearch(@(par) negL(par,exceed,U(k)),init);
    init_temp = init;
    while_counter = 0;
    while param == init_temp                                                    % in case initial guess is bad
        while_counter = while_counter + 1;
        init_temp = [max(0.1, init(1) + normrnd(0,shake_guess^2)), init(2) + normrnd(0,shake_guess^2)];
        param = fminsearch(@(par) negL(par,exceed,U(k)),init_temp);
        if while_counter == 10
            break
            'pick better initial value'
        end
    end
    param_save(:,k) = param;
    p_u = length(exceed)/length(data);
    pc(k) = p_u*(1 - gpcdf(trans(0), param(2), param(1),U(k)) );
    ue = U(k) - param(1)/param(2);
    if param(2)<0; ue_save(k) = ue; end

    if qqplot == 1

        %%% generate sample of stochastic ttc, where one sample is drawn from each encounter.
        col_ind = randsample(Nexp,Nenc, true)';
        ind = sub2ind(size(trans_data), 1:Nenc, col_ind );
        ttc_sample = trans_data(ind);
        exceed = ttc_sample(ttc_sample>U(k));
        excess = exceed - U(k);

        %%% compute empirical d.f. for excesses
        exceed_data = trans_data(   max(trans_data,[],2) > U(k) ,:   );   % keeps only rows with atleast one obs. above threshold
        pu_i = sum(exceed_data>U(k), 2)/ size(trans_data,2);                     % probability to exceed threshold for each encounter
        weights = pu_i/sum(pu_i) ;

        % replace obs. below u with nan's
        below_u = find(exceed_data<=U(k));
        exceed_data(below_u) = nan*ones(1,length(below_u));
        exceed_data = exceed_data - U(k); % subtract u to get excessses

        % evaluate empirical distribution function
        xplot_lower = 0;
        xplot_upper = max(excess)*1.1;
        x_eval = linspace(xplot_lower, xplot_upper, n_eval_cdf);
        femp = weights'*F_emp(x_eval, exceed_data);



        clf
        subplot(211)
            qq_plot(exceed,param(1),param(2),U(k),k)
        subplot(212)
            plot(x_eval,femp,'.')
            hold on
            plot(x_eval, gpcdf(x_eval, param(2), param(1), 0))
            line(trans([0,0]), 1.2,'LineStyle','--');
            line(get(gca, 'xlim'), [1 1],'Color','green','LineStyle','--');
            title(sprintf('empirical d.f. vs model d.f. for threshold %d.',k))

        pause(pause_length)
    end

    init = param;

%%%%%%%%%%%%%%%%%%%%%%%%%% bootstrapping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if compute_ci == 1
        sigma_sample = zeros(1,Nbs);
        xi_sample = zeros(1,Nbs);
        p_nea_sample = zeros(1,Nbs);
        for j=1:Nbs

            resampling = randsample(Nenc,Nenc,true);              % indeces used to bootstrap
            trans_DAFEA_bs = trans_data(resampling,:);            % bootstrapped sample
            data = trans_DAFEA_bs(:);
            exceed = data(data>U(k));                         % get exceedences
            param_bs =  fminsearch(@(par) negL(par,exceed,U(k)),init);                           % estimate parameters
            init_temp = init;
            while_counter = 1;

            while param_bs == init_temp
                while_counter = while_counter + 1;
                % in case initial guess is bad
                init_temp = [max(0.1, init(1) + normrnd(0,shake_guess^2)), init(2) + normrnd(0,shake_guess^2)];
                param_bs = fminsearch(@(par) negL(par,exceed,U(k)),init_temp);

                if while_counter==300; param_bs = [nan,nan]; break; end
            end
            p_u = length(exceed)/length(data);
            sigma_sample(j) = param_bs(1); %#ok<*SAGROW>
            xi_sample(j) = param_bs(2);
            p_nea_sample(j) = p_u*(1 - gpcdf(trans(0), param_bs(2), param_bs(1),U(k)) );

        end
        sigma_sample(isnan(sigma_sample)) = [];
        xi_sample(isnan(xi_sample)) = [];
        p_nea_sample(isnan(p_nea_sample)) = [];

        sigma_mean = mean(sigma_sample);
        xi_mean = mean(xi_sample);
        p_nea_mean = mean(p_nea_sample);

        se_sigma = sqrt(sum((sigma_sample - sigma_mean).^2)/length(sigma_sample));
        se_xi = sqrt(sum((xi_sample - xi_mean).^2)/length(xi_sample));
        se_p_nea = sqrt(sum((p_nea_sample - p_nea_mean).^2)/length(p_nea_sample));

        ci_sigma = param(1) + [-1 1]*1.96*se_sigma; ci_sigma = sort(ci_sigma);
        ci_xi = param(2) + [-1 1]*1.96*se_xi; ci_xi = sort(ci_xi);
        ci_p_nea = pc(k) + [-1 1]*1.96*se_p_nea; ci_p_nea = sort(ci_p_nea);

        ci_sigma_u(:,k) = ci_sigma';
        ci_xi_u(:,k) = ci_xi';
        ci_p_nea_u(:,k) = ci_p_nea';

    end

end

clf;
subplot(221)
plot(U,param_save(1,:))
title('sigma_{est}')
hold on
if compute_ci == 1
    plot(U,ci_sigma_u)
end
subplot(222)
plot(U, param_save(2,:),'s'); hold on
plot(U, param_save(2,:),'b')
xlabel('threshold')
title('xi_{est}')
if compute_ci == 1
    plot(U,ci_xi_u,':','color','b')
end

subplot(223)
if min(param_save(2,:))<0
    plot(U,ue_save); hold on
    plot(U,ones(1,m)*max_data)
    plot(U,ones(1,m)*trans(0))
    title('upper endpoint estimates')
end
subplot(224)
% if logit==1
%     plot(U,log10(pc))
%     title('log(p_{est})')
%     if compute_ci == 1
%         hold on
%         %plot(U,ci_p_nea_u)
%     end
% else
%     plot(U,pc)
%     title('p_{est}')
%     if compute_ci == 1
%         hold on
%         %plot(U,ci_p_nea_u)
%     end
% end
thr_autofind(ci_xi_u, param_save(2,:),m)

%% estimating p(collision type i)
p_interactive = (sum(enc_type==-1) + sum(enc_type==-2) + sum(enc_type==2))/N; % probability of encounter being interactive
p_ea = (sum(enc_type==-1)+sum(enc_type==-2))/N
% this estimator gives P(C1,NEA) when p_nea is based on danger at first
% interactive action, DAFEA or danger_FEA
if data_type == 1
    p_c1 = p_interactive*pc
end

% this estimaor gives P(C, NEA) when p_nea is based on X
if data_type == 2
    p_c1 = pc
end

% this estimaor gives P(C, EA) when p_nea is based on X
if data_type == 3
    p_c2 = pc*p_ea
end
