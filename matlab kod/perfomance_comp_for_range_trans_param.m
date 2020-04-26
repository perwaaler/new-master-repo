
%% Estimation of P(collision) using different transformation parameters
% this code is used to take simulated data and perform POT method to
% estimate P(collision). Data_matrix  can be either danger_FEA, DAFEA, X,
% or danger_max_EA. data_type=1 --> danger_FEA or DAFEA, data_type=2 --> X,
% data_type=3 --> danger_max_EA
% 1=DAFEA 2=X 3=ttc_FEA 4=ttc_min ttc_min_EA=5 danger_FEA=6 dist_min_EA=7 dist_min=8

data_type_list=[2]; %#ok<NBRAK>
select_trans = 3;
% par_range = [.1 .3 .5 .7 .9 1.2 1.4 1.5];
par_range = [5.5];
% par_range = [5];
% par_range = [.1 .2 .3 .4 .5];
if select_trans==1; par_range = 1;end
safety_level = 1;

% maximum and minimum fraction of data to use after applying thresholds
up_frac = 0.8;
lo_frac = 0.06;

% plotting options
compute_ci = 0;              % set equal to one if confidence intervals for xi are desired
qqplot = 0;
save_plot = 0;

% pausing options
pause_trans = 0.0;
qq_pause = 0.0;
stability_pause = 0;

for mm=1:length(data_type_list)
    
data_type = data_type_list(mm); % 1=DAFEA 2=X 3=ttc_FEA 4=ttc_min ttc_min_EA=5 danger_FEA=6 dist_min_EA=7 dist_min=8

for i=1:length(par_range)
    % transformation parameters
    if select_trans==2
        p_ex = 0.1;%[0.01 .03 .05 .07 .09];
        p_ex = p_ex(i);
        trans_par = p_ex;
    else
        p_inv = par_range(i);%3;%[0.1 .3 .5 .7 .9]; %1 + (i - 1)/2;
        d_inv = 3.5;
        trans_par = [p_inv, d_inv];
    end

    if safety_level==1
        load data_500enc_r_0_3_safety_level_1.mat
    else
        load data_500enc_r_0_3_safety_level_2.mat
    end

    J = size(all_data,1);

    % number of thresholds
    n_u = 10;

    % arrays for saving data
    ue_save_matrix = zeros(J,n_u);
    param_save_matrix = zeros(J,n_u);
    pc_save_matrix = zeros(J,n_u);
    p_c_save_matrix = zeros(J,n_u);
    ci_xi_u_matrix = cell(J, 1);
    thr_save_matrix = zeros(J,n_u);
    p_exceed_matrix = zeros(J,n_u);

    % inital value for parameter estimation
    init_val = [1 0];
    
    trans =@(x) transform(x, select_trans, trans_par);

    for jj=1:J
        % generate sample of encounters
        jj %#ok<NOPTS>

        data_matrix = all_data{jj,data_type}; % select data. row i should correspond to encounter i, and column j to j'th simulated ttc value (in case of stochastic ttc)

        % find encounters with finite ttc values
        min_data = data_matrix(min(data_matrix,[],2)<Inf,:);
        min_data = min(min_data,[],2);

        % find minimum and maximum thresholds based on amount of data to be used
        u_minmax = find_threshold(min_data, lo_frac, up_frac);
        u_l = u_minmax(1);
        u_u = u_minmax(2);

        % transform thresholds and data
        u_min_max = sort(trans([u_l,u_u]));
        u_l_trans = u_min_max(1);
        u_u_trans = u_min_max(2);
        U = sort(trans(linspace(u_l, u_u, n_u)));
        thr_save_matrix(jj,:) = U;
        trans_data = trans(data_matrix);

        % plot transformed data with transformed thresholds
%         clf
%         plot(trans(min_data),'.')
%         hold on
%         plot(1:length(min_data),U'*ones(1,length(min_data)), 'k')
%         plot(ones(1,length(min_data))*trans(0))
%         ylim([trans(30),trans(0)*1.1])
%         title('transformed data')
%         pause(pause_trans)

        %%% ensure good initial value!
        init = est_par(trans_data(:), init_val, U(1));

        %%% Fit models for range of thresholds and show diagnostic plots
        shake_var = 0.1;                 % variance of noise that gets added to initial guess when stuck
        Nenc = length(data_matrix(:,1));   % number of encounters
        Nexp = length(data_matrix(1,:));   % number of values per row. If surrogate-measure is deterministic, then Nexp=1.
        Nbs = 20;                         % number of bootstrapped samples to compute standard error

        % collects 95% ci's for xi and sigma
        ci_sigma_u = zeros(2,n_u)*nan;
        ci_xi_u = zeros(2,n_u)*nan;                                      % collects 95% ci's for xi
        ci_p_nea_u = zeros(2,n_u)*nan;

        param_save = zeros(2,n_u);
        pc = zeros(1,n_u)*nan;                                                           % collects estimated collision probability for each threshold
        ue_save = zeros(1,n_u)*nan;                                                      % collects estimated upper endpoint
        max_data = max(max(trans_data));                                                 % largest observed value
        logit = 1;

        p_inter = all_data{jj,end-1};
        p_ea = all_data{jj,end};

        % limits for plots of empirical and model distribution functions
        xplot_lower = 0;
        xplot_upper = 1;

        % number of points where cdf gets evaluated
        n_eval_cdf = 1000;

        for k=1:n_u

            param = est_par_w_shake(trans_data(:), init, U(k), shake_var);
            param_save(:,k) = param;
            p_u = sum(trans_data(:) > U(k))/(p_inter*500*size(trans_data,2));
            p_exceed_matrix(jj,k) = p_u;

            pc(k) = p_u*(1 - gpcdf(trans(0), param(2), param(1), U(k)) );
            ue = U(k) - param(1)/param(2);
            if param(2)<0; ue_save(k) = ue; end

            if qqplot == 1

                %%% generate sample of stochastic ttc, where one sample is drawn from each encounter.
                col_ind = randsample(Nexp, Nenc, true)';
                ind = sub2ind(size(trans_data), 1:Nenc, col_ind );
                ttc_sample = trans_data(ind);
                exceed = ttc_sample( ttc_sample>U(k) );
                excess = exceed - U(k);

                %%% compute empirical d.f. for excesses
                exceed_data = trans_data(    max(trans_data,[],2) > U(k) ,:   );   % keeps only rows with atleast one obs. above threshold
                pu_i = sum(exceed_data>U(k), 2)/ size(trans_data,2);                     % probability to exceed threshold for each encounter
                weights = pu_i/sum(pu_i) ;

                %%% replace obs. below u with nan's
                below_u = find(exceed_data<=U(k));
                exceed_data(below_u) = nan*ones(1,length(below_u));
                exceed_data = exceed_data - U(k); % subtract u to get excessses

                %%% evaluate empirical distribution function
                xplot_lower = 0;
                xplot_upper = max(excess)*3;
                x_eval = linspace(xplot_lower, xplot_upper, n_eval_cdf);
                femp = weights'*F_emp(x_eval, exceed_data);

                %%% plot distribution functions and qq-plots
                clf
                subplot(211)
                    qq_plot(exceed,param(1),param(2),U(k),k)
                    title(sprintf('empirical vs model quantiles, for threshold %d.',k))
                subplot(212)
                    plot(x_eval,femp,'.')
                    hold on
                    plot(x_eval, gpcdf(x_eval, param(2), param(1), 0))
                    line(([0,0]), 1.2,'LineStyle','--');
                    line(get(gca, 'xlim'), [1 1],'Color','green','LineStyle','--');
                    title(sprintf('empirical d.f. vs model d.f. for threshold %d.',k))
                % save plot
                if save_plot == 1
                    savefig(sprintf('goodness_of_fit_test_DAFEA_exp_%d_%d_id4',jj,k))
                end
                pause(qq_pause)
            end

            init = param;

            %%%%%%%%%%%%%%%%%% bootstrapping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if compute_ci == 1

                sigma_sample = zeros(1,Nbs);
                xi_sample = zeros(1,Nbs);
                p_nea_sample = zeros(1,Nbs);

                for j=1:Nbs

                    [param_bs, p_u_bs] = bootstrap_est(trans_data, Nenc, init, U(k), shake_var);
                    sigma_sample(j) = param_bs(1); %#ok<*SAGROW>
                    xi_sample(j) = param_bs(2);
                    p_nea_sample(j) = p_u_bs*(1 - gpcdf(trans(0), param_bs(2), param_bs(1),U(k)) );

                end

                [ci_sigma, ci_xi, ci_pc] = compute_bs_ci(param, pc(k), sigma_sample, xi_sample, p_nea_sample);
                ci_sigma_u(:,k) = ci_sigma';
                ci_xi_u(:,k) = ci_xi';
                ci_p_nea_u(:,k) = ci_pc';


            end

        end

        %%%%% diagnostic plots %%%%%
        clf;
        subplot(221)%%%%%%%%%%%%%%%%%%%%%%
        plot(U,param_save(1,:))
        title('sigma_{est}')
        hold on
        if compute_ci == 1
            plot(U,ci_sigma_u)
        end

        subplot(222)%%%%%%%%%%%%%%%%%%%%%%
        plot(U   , param_save(2,:))
        title('xi_{est}')
        hold on
        if compute_ci == 1
            plot(U,ci_xi_u)
            thr_index = thr_autofind(ci_xi_u, param_save(2,:), n_u);
            plot(U(thr_index), param_save(2,thr_index),'*');
        end

        subplot(223)%%%%%%%%%%%%%%%%%%%%%%
        if min(param_save(2,:))<0
            plot(U,ue_save); hold on
            plot(U,ones(1,n_u)*max_data)
            plot(U,ones(1,n_u)*trans(0))
            title('upper endpoint estimates')
        end
        subplot(224)%%%%%%%%%%%%%%%%%%%%%%
        if logit==1
            plot(U,log10(pc))
            title('log(p_{est})')
            if compute_ci == 1
                hold on
                %plot(U,ci_p_nea_u)
            end
        else
            plot(U,pc)
            title('p_{est}')
            if compute_ci == 1
                hold on
                %plot(U,ci_p_nea_u)
            end
        end

        param_save_matrix(jj,:) = param_save(1,:) + 1i*param_save(2,:);
        ue_save_matrix(jj,:) = ue_save;
        pc_save_matrix(jj,:) = pc;
        ci_xi_u_matrix{jj,1} = ci_xi_u;

        pause(stability_pause)

        %%% estimating p(collision) %%%

        % this estimator gives P(C,NEA) when p_nea is based on severity
        % measure at first evasive action.
        if data_type == 1 || data_type == 3 || data_type == 6
            p_c = p_inter*pc;
        end

        % this estimaor gives P(C, NEA) when p_nea is based on maximum
        % or minimum severity measure over full encounter
        if data_type == 2 || data_type == 4 || data_type == 8
            p_c = pc;
        end

        % this estimaor gives P(C, EA) when p_nea is based on maximum
        % or minimum severity level over attempt to avoid collision
        if data_type == 5 || data_type == 7
            p_c = pc*p_ea;
        end

        p_c_save_matrix(jj,:) = p_c;


    end

    % save information
    hit_rate = sum(pc_save_matrix>0,1)/500*100;

    % save results
    variables = cell(1,7);
    variables{1,1} = hit_rate;
    variables{1,2} = pc_save_matrix;
    variables{1,3} = p_c_save_matrix;
    variables{1,4} = thr_save_matrix;
    variables{1,5} = param_save_matrix;
    variables{1,6} = ci_xi_u_matrix;
    variables{1,7} = p_exceed_matrix;
    
    save_result(variables, data_type, select_trans, trans_par, safety_level, up_frac, lo_frac)
    
end






end
