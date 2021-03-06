%% generate data
E = 6;
V = 10;




% parametes as fcn of mean and variance
a = E^2/V;
b = V/E;
alfa = 0.078;
a = 0.5;
b = 1;

normcdf(0,3.5,.9)
chi2cdf(alfa,9)


data_Matrix = 1./gamrnd(a*ones(500,500),b*ones(500,500))-alfa;
% data_Matrix = chi2rnd(9*ones(500,500))-alfa;
% data_Matrix = betarnd(ones(500,500)*5,ones(500,500)*3)-alfa;
data_Matrix = normrnd(3.5*ones(500,500),.9*ones(500,500));
p_actual = gamcdf(alfa,a,b);
% p_actual = chi2cdf(alfa,9);
% p_actual = betacdf(alfa,5,3)
p_actual = normcdf(0,3.5,.9)

clf
plot(data_Matrix(1,:),'.')
ylim([-1,20])
%%
%plot(linspace(0,20,100), normpdf(linspace(0,20,100),3.5,1))
plot(linspace(0,20,100), chi2pdf(linspace(0 ,20,100),15))
plot(linspace(0,1,100), betapdf(linspace(0 ,1,100),5,3))
plot(linspace(-1,3,1000), normpdf(linspace(-1 ,3,1000),3.5,.9))
%%
plot_diagnostics = 0;
pEst_container = cell(7,1);
pEst_matrix = zeros(500,10);

for J=4:7
J
pEst_matrix = zeros(500,10);
compute_ci = 0;                    % set equal to one if confidence intervals for xi are desired
qqplot = 0;
qq_pause = .5;

for jj = 1:500
    jj; %#ok<NOPTS>

    data_matrix = data_Matrix(:,jj);
    up_frac = 0.50;
    lo_frac = 0.04;
    safety_level = 1;
    data_type = 2;
    sev_ind = sev_measure(data_type);

    % find encounters with finite ttc values
    min_data = data_matrix(min(data_matrix,[],2)<Inf,:);
    min_data = min(min_data,[],2);
    save_plot = 0;% find minimum and maximum thresholds based on amount of data to be used

    u_minmax = find_threshold(min_data, lo_frac, up_frac);

    u_l = u_minmax(1);
    u_u = u_minmax(2);

%     clf;
%     plot(min_data,'.'); hold on 
%     plot(ones(1,length(min_data))*u_l) 
%     plot(ones(1,length(min_data))*u_u)
%     %ylim([-1,100])
%     title('untransformed data')
%     transforming data and plotting transformed data and thresholds

    select_trans = 2;    

    % transformation choice and parameters
    if select_trans==1
        trans_par = [];
    elseif select_trans==2
        trans_par = 0.4*J;%select_par(kk);%p_ex;
    else
        trans_par = [2.5, 3.5];
    end

<<<<<<< HEAD
up_frac = 0.80;
lo_frac = 0.06;
safety_level=1;
data_type = 2;
sev_ind = sev_measure(data_type);

% find encounters with finite ttc values
min_data = data_matrix(min(data_matrix,[],2)<Inf,:);
min_data = min(min_data,[],2);
save_plot = 0;% find minimum and maximum thresholds based on amount of data to be used

u_minmax = find_threshold(min_data, lo_frac, up_frac);

u_l = u_minmax(1);
u_u = u_minmax(2);

% clf;
% plot(min_data,'.'); hold on 
% plot(ones(1,length(min_data))*u_l) 
% plot(ones(1,length(min_data))*u_u)
% %ylim([-1,100])
% title('untransformed data')
% transforming data and plotting transformed data and thresholds

select_trans =1;    
 
% transformation choice and parameters
if select_trans==1
    trans_par = [];
elseif select_trans==2
    trans_par = 0.30;%select_par(kk);%p_ex;
else
    trans_par = [2.5, 3.5];
end

trans =@(x) transform(x, select_trans, trans_par);

% transform thresholds
u_min_max = sort(trans([u_l,u_u]));
u_l_trans = u_min_max(1);
u_u_trans = u_min_max(2);
m=10;
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
init = fminsearch(@(par) negL(par, exceed, u_l_trans), [3 0.3]);

% Fit models for range of thresholds and show diagnostic plots
shake_guess = 0.1;                 % variance of noise that gets added to initial guess when stuck
Nenc = length(data_matrix(:,1));   % number of encounters
Nexp = length(data_matrix(1,:));   % number of values per row. If surrogate-measure is deterministic, then Nexp=1.
Nbs = 200;                          % number of bootstrapped samples to compute standard error
trans_data = trans(data_matrix);                                                  % number of thresholds used for estimation
                   % vector containing thresholds'

% collects 95% ci's for xi
ci_sigma_u = zeros(2,m)*nan;
ci_xi_u = zeros(2,m)*nan;                                      % collects 95% ci's for xi
ci_p_nea_u = zeros(2,m)*nan;
se_p_nea_save = zeros(1,m);

param_save = zeros(2,m);
pc = zeros(1,m)*nan;                                    % collects estimated collision probability for each threshold
ue_save = zeros(1,m)*nan;                                  % collects estimated upper endpoint
max_data = max(max(trans_data));                              % largest observed value
negL = @(par, exceed_data,u) -sum( log(gppdf(exceed_data,par(2),par(1),u)) );  %negative log likelihood fcn.
logit = 1;                                               % set to plot logarithm of p_nea when magnitude of p_nea varies alot
compute_ci = 1;                    % set equal to one if confidence intervals for xi are desired
qqplot = 1;
qq_pause = 0;
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
    ue_save(k) = U(k) - param(1)/param(2);
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
        xplot_upper = max(excess*1.4);
        x_eval = linspace(xplot_lower, xplot_upper, n_eval_cdf);
        femp = weights'*F_emp(x_eval, exceed_data);



        
%        sgtitle(sprintf('goodness of fit plots, %s %s %s, threshold %d',sevme_str(sev_ind), trans_str(select_trans),num2str(trans_par), k))
        clf        
        subplot(211)
            qq_plot(exceed,param(1),param(2),U(k),k)
            title('empirical vs model quantiles')
        subplot(212)
            
            plot(x_eval,femp,'.')
            hold on
            plot(x_eval, gpcdf(x_eval, param(2), param(1), 0))
            %line(trans([0,0]), 1.2,'LineStyle','--');
            line(get(gca, 'xlim'), [1 1],'Color','green','LineStyle','--');
            %title('empirical vs model distribution functions')
        if save_plot==1
            saveas(gcf, sprintf('goodness_of_fit_sample_%d_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d_trhind_%d.png',sample, data_type,select_trans,100*trans_par,safety_level,up_frac*100,lo_frac*100,k))
            %savefig(sprintf('goodness_of_fit_mindistFEA_thr_%d',k))
=======
    trans =@(x) transform(x, select_trans, trans_par);

    % transform thresholds
    u_min_max = sort(trans([u_l,u_u]));
    u_l_trans = u_min_max(1);
    u_u_trans = u_min_max(2);
    m = 10;
    U = sort(trans(linspace(u_l, u_u,m)));

    trans_data = trans(data_matrix);

    clf; plot(trans(min_data),'.'); hold on
    plot(1:length(min_data),U'*ones(1,length(min_data)), 'k');plot(ones(1,length(min_data))*trans(0));ylim([trans(30),trans(0)*1.1])
    title('transformed data')
    legend('transformed thresholds')
    % ensure good initial value!
    trans_data = trans_data(:);
    exceed = trans_data(trans_data > u_l_trans);
    negL = @(par,exceed_data,u) -sum( log(gppdf(exceed_data,par(2),par(1),u)) );
    init = fminsearch(@(par) negL(par, exceed, u_l_trans), [3 0.3]);

    % Fit models for range of thresholds and show diagnostic plots
    shake_guess = 0.1;                 % variance of noise that gets added to initial guess when stuck
    Nenc = length(data_matrix(:,1));   % number of encounters
    Nexp = length(data_matrix(1,:));   % number of values per row. If surrogate-measure is deterministic, then Nexp=1.
    Nbs = 200;                          % number of bootstrapped samples to compute standard error
    trans_data = trans(data_matrix);                                                  % number of thresholds used for estimation
                       % vector containing thresholds'

    % collects 95% ci's for xi
    ci_sigma_u = zeros(2,m)*nan;
    ci_xi_u = zeros(2,m)*nan;                                      % collects 95% ci's for xi
    ci_p_nea_u = zeros(2,m)*nan;
    se_p_nea_save = zeros(1,m);

    param_save = zeros(2,m);
    pc = zeros(1,m)*nan;                                    % collects estimated collision probability for each threshold
    ue_save = zeros(1,m)*nan;                                  % collects estimated upper endpoint
    max_data = max(max(trans_data));                              % largest observed value
    negL = @(par, exceed_data,u) -sum( log(gppdf(exceed_data,par(2),par(1),u)) );  %negative log likelihood fcn.
    logit = 1;                                               % set to plot logarithm of p_nea when magnitude of p_nea varies alot

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
            xplot_upper = max(excess*1.01);
            x_eval = linspace(xplot_lower, xplot_upper, n_eval_cdf);
            femp = weights'*F_emp(x_eval, exceed_data);




    %        sgtitle(sprintf('goodness of fit plots, %s %s %s, threshold %d',sevme_str(sev_ind), trans_str(select_trans),num2str(trans_par), k))
    %         subplot(211)
    %             qq_plot(exceed,param(1),param(2),U(k),k)
    %             title('empirical vs model quantiles')

                clf
                plot(x_eval,femp,'.')
                hold on
                plot(x_eval, gpcdf(x_eval, param(2), param(1), 0))
                %line(trans([0,0]), 1.2,'LineStyle','--');
                line(get(gca, 'xlim'), [1 1],'Color','green','LineStyle','--');
                %title('empirical vs model distribution functions')
            if save_plot==1
                saveas(gcf, sprintf('goodness_of_fit_sample_%d_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d_trhind_%d.png',sample, data_type,select_trans,100*trans_par,safety_level,up_frac*100,lo_frac*100,k))
                %savefig(sprintf('goodness_of_fit_mindistFEA_thr_%d',k))
            end

            pause(qq_pause)
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
            se_p_nea_save(k) = se_p_nea(1);

            ci_sigma = param(1) + [-1 1]*1.96*se_sigma; ci_sigma = sort(ci_sigma);
            ci_xi = param(2) + [-1 1]*1.96*se_xi; ci_xi = sort(ci_xi);
            ci_p_nea = pc(k) + [-1 1]*1.96*se_p_nea; ci_p_nea = sort(ci_p_nea);

            ci_sigma_u(:,k) = ci_sigma';
            ci_xi_u(:,k) = ci_xi';
            ci_p_nea_u(:,k) = ci_p_nea';

>>>>>>> 12249eb4234c713dfd193ae48745062510f013e6
        end

    end
    %standard_error{sample,kk} = se_p_nea_save;



    clf
    %sgtitle(sprintf('shape paremter and probability estimates, %s %s %s',sevme_str(sev_ind), trans_str(select_trans),num2str(trans_par)))
    subplot(211)
    plot(U, param_save(2,:),'s')
    hold on
    plot(U, param_save(2,:),'b')
    xlabel('threshold')
    title('xi_{est}')
    if compute_ci == 1
        plot(U,ci_xi_u,':','color','b')
    end

    if plot_diagnostics==1
        subplot(212)
        if logit==1
            plot(U,log10(pc),'s')
            hold on
            plot(U,log10(pc))
            title('log_{10}(p_{est})')
            xlabel('threshold')
            if compute_ci == 1
                hold on
                %plot(U,ci_p_nea_u)
            end
        else
            plot(U,pc,'s')
            hold on
            plot(U,pc)
            title('p_{est}')
            if compute_ci == 1
                hold on
                plot(U,ci_p_nea_u)
            end

        end
    end
    guess(select_trans)=10^mean(log10(pc));
    pEst_matrix(jj,:) = pc;
    end
pEst_container{J,1} = pEst_matrix;
if select_trans == 1
    break
end
%%%%%%%%%%%%%%% plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% plot results for gamma distribution
%%

SR_thr = 1;
SR_NoTrans = mean(abs(log10(pEst_matrix_noTrans/p_actual))<SR_thr)*100;
SR_10 = mean(abs(log10(pEst_matrix_10/p_actual))<SR_thr)*100;
SR_20 = mean(abs(log10(pEst_matrix_20/p_actual))<SR_thr)*100;
SR_30 = mean(abs(log10(pEst_matrix_30/p_actual))<SR_thr)*100;
SR_40 = mean(abs(log10(pEst_matrix_40/p_actual))<SR_thr)*100;
SR_50 = mean(abs(log10(pEst_matrix_50/p_actual))<SR_thr)*100;
SR_60 = mean(abs(log10(pEst_matrix_60/p_actual))<SR_thr)*100;
SR_70 = mean(abs(log10(pEst_matrix_70/p_actual))<SR_thr)*100;

clf
color = linspace(0,1,7);
plot(SR_NoTrans,'red')
hold on
plot(SR_10, 'marker','o', 'color', [0 color(1) 1])
plot(SR_20, 'marker', '*','color', [0 color(2) 1])
plot(SR_30,'marker', '+','color', [0 color(3) 1])
plot(SR_40, 'marker', 'x','color', [0 color(4) 1])
plot(SR_50, 'marker', 'p','color', [0 color(5) 1])
plot(SR_60, 'marker', 'd','color', [0 color(6) 1])
plot(SR_70, 'marker', 'h','color', [0 color(7) 1])

legend('No Trans','p = 0.1', 'p = 0.2','p = 0.3','p = 0.4','p = 0.5','p = 0.6','p = 0.7')

% exp(0.4) 10^-4

%% plot results for chisquared distribution



SR_thr = 1;
SR_NoTrans = mean(abs(log10(pEst_matrix_noTrans/p_actual))<SR_thr)*100;
SR_10 = mean(abs(log10(pEst_container{1,1}/p_actual))<SR_thr)*100;
SR_20 = mean(abs(log10(pEst_container{2,1}/p_actual))<SR_thr)*100;
SR_30 = mean(abs(log10(pEst_container{3,1}/p_actual))<SR_thr)*100;
SR_40 = mean(abs(log10(pEst_container{4,1}/p_actual))<SR_thr)*100;
SR_50 = mean(abs(log10(pEst_container{5,1}/p_actual))<SR_thr)*100;
SR_60 = mean(abs(log10(pEst_container{6,1}/p_actual))<SR_thr)*100;
SR_70 = mean(abs(log10(pEst_container{7,1}/p_actual))<SR_thr)*100;
% SR_80 = mean(abs(log10(pEst_container{8,1}/p_actual))<SR_thr)*100;
% SR_100 = mean(abs(log10(pEst_container{10,1}/p_actual))<SR_thr)*100;
% SR_110 = mean(abs(log10(pEst_container{11,1}/p_actual))<SR_thr)*100;
clf
color = linspace(0,1,7);
plot(SR_NoTrans,'red')
hold on
plot(SR_10, 'marker','o', 'color', [0 color(1) 1])
plot(SR_20, 'marker', '*','color', [0 color(2) 1])
plot(SR_30,'marker', '+','color', [0 color(3) 1])
plot(SR_40, 'marker', 'x','color', [0 color(4) 1])
plot(SR_50, 'marker', 'p','color', [0 color(5) 1])
plot(SR_60, 'marker', 'd','color', [0 color(6) 1])
plot(SR_70, 'marker', 'h','color', [0 color(7) 1])
% plot(SR_80, 'marker', '^','color', [0 color(8) 1])
% plot(SR_100, 'marker', '>','color', [0 color(10) 1])
% plot(SR_110, 'marker', '<','color', [0 color(11) 1])

legend('No Trans','p = 0.5', 'p = 1.0','p = 1.5','p = 2.0','p = 2.5','p = 3.0','p = 3.5','p = 4.0','4.5')

%% plot results for normal distribution



SR_thr = 1;
SR_NoTrans = mean(abs(log10(pEst_matrix_noTrans/p_actual))<SR_thr)*100;
SR_10 = mean(abs(log10(pEst_container{1,1}/p_actual))<SR_thr)*100;
SR_20 = mean(abs(log10(pEst_container{2,1}/p_actual))<SR_thr)*100;
SR_30 = mean(abs(log10(pEst_container{3,1}/p_actual))<SR_thr)*100;
SR_40 = mean(abs(log10(pEst_container{4,1}/p_actual))<SR_thr)*100;
SR_50 = mean(abs(log10(pEst_container{5,1}/p_actual))<SR_thr)*100;
SR_60 = mean(abs(log10(pEst_container{6,1}/p_actual))<SR_thr)*100;
SR_70 = mean(abs(log10(pEst_container{7,1}/p_actual))<SR_thr)*100;
% SR_80 = mean(abs(log10(pEst_container{8,1}/p_actual))<SR_thr)*100;
% SR_100 = mean(abs(log10(pEst_container{10,1}/p_actual))<SR_thr)*100;
% SR_110 = mean(abs(log10(pEst_container{11,1}/p_actual))<SR_thr)*100;
clf
color = linspace(0,1,7);
plot(SR_NoTrans,'red')
hold on
plot(SR_10, 'marker','o', 'color', [0 color(1) 1])
plot(SR_20, 'marker', '*','color', [0 color(2) 1])
plot(SR_30,'marker', '+','color', [0 color(3) 1])
plot(SR_40, 'marker', 'x','color', [0 color(4) 1])
plot(SR_50, 'marker', 'p','color', [0 color(5) 1])
plot(SR_60, 'marker', 'd','color', [0 color(6) 1])
plot(SR_70, 'marker', 'h','color', [0 color(7) 1])
% plot(SR_80, 'marker', '^','color', [0 color(8) 1])
% plot(SR_100, 'marker', '>','color', [0 color(10) 1])
% plot(SR_110, 'marker', '<','color', [0 color(11) 1])

legend('No Trans','p = 0.4', 'p = 0.8','p = 1.2','p = 1.6','p = 2.0','p = 2.4','p = 2.8')






