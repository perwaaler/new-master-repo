%% Estimation of P(collision) using different transformation parameters
% this code is used to take simulated data and perform POT method to
% estimate P(collision). Data_matrix  can be either danger_FEA, DAFEA, X,
% or danger_max_EA. data_type=1 --> danger_FEA or DAFEA, data_type=2 --> X,
% data_type=3 --> danger_max_EA

for safety_level=1:2
    
    % select which severity measure to use
    data_type = 1;          
    
    % select data set: 2 --> comparison set, 1 --> default settings
    if safety_level == 1
        load data_500enc_r_safety_level_1.mat
    else
        load data_500enc_r_safety_level_2.mat
    end

    %%%%% loop over different transformation parameters %%%%%
    for ii=(1:5)*2
        ii %#ok<NOPTS>
        n_u = 10; 

        % number of encounter samples
        J = size(all_data,1);              
                                                                   % number of thresholds used for estimation
        % what fraction of the data to use
        upper_fraction = 0.8;
        lower_fraction = 0.15;
        
        % inital value for parameter estimation
        init_val = [3 0.3];

        % transformation selection and parameters
        select_trans = 2;
        p_ex = ii/10;
        p_inv = 3;
        d_inv = 3;
        trans =@(x) transform(x, select_trans, p_ex);

        % plotting options
        compute_ci = 1;                                       % set equal to one if confidence intervals for xi are desired
        qqplot = 1;
        save_plot = 0;

        % arrays for saving data                                                         
        ue_save_matrix = zeros(J,n_u);
        param_save_matrix = zeros(J,n_u);
        pc_save_matrix = zeros(J,n_u);
        p_c_save_matrix = zeros(J,n_u);
        ci_xi_u_matrix = cell(J, 1); 
        thr_save_matrix = zeros(J,n_u);
        p_exceed_matrix = zeros(J,n_u);

        % pausing options
        pause_trans = 0.0;
        qq_pause = 0.0;
        stability_pause = 0.0;
        
            %%%%% loop over each all sets of encounters %%%%%
        for jj=1:J
            jj %#ok<NOPTS>
            p_inter = all_data{jj,end-1};
            p_ea = all_data{jj,end};

            % generate sample of encounters

            data_matrix = all_data{jj,data_type}; % select data. row i should correspond to encounter i, and column j to j'th simulated ttc value (in case of stochastic ttc)

            % find encounters with finite ttc values
            min_data = data_matrix( min(data_matrix,[],2)<Inf ,:);
            min_data = min(min_data,[],2);

            % find minimum and maximum thresholds based on amount of data to be used
            u_minmax = find_threshold(min_data, lower_fraction, upper_fraction);
            u_l = u_minmax(1);
            u_u = u_minmax(2);

            clf;plot(min_data,'.'); hold on; plot(ones(1,length(min_data))*u_l); plot(ones(1,length(min_data))*u_u)
            ylim([-1,30])

            % transform thresholds and data
            u_min_max = sort(trans([u_l,u_u]));
            u_l_trans = u_min_max(1); 
            u_u_trans = u_min_max(2);
            U = sort(trans(linspace(u_l, u_u, n_u)));
            thr_save_matrix(jj,:) = U;
            trans_data = trans(data_matrix);
            
            % plot transformed data and thresholds
            clf; plot(trans(min_data),'.'); hold on
            plot(1:length(min_data),U'*ones(1,length(min_data)), 'k');plot(ones(1,length(min_data))*trans(0));ylim([trans(30),trans(0)*1.1])
            title('transformed data')
            pause(pause_trans)

            %%% ensure good initial value!
            init = est_par(trans_data(:), init_val, U(1));

            %%% Fit models for range of thresholds and show diagnostic plots
            shake_guess = 0.1;                 % variance of noise that gets added to initial guess when stuck
            Nenc = length(data_matrix(:,1));   % number of encounters
            Nexp = length(data_matrix(1,:));   % number of values per row. If surrogate-measure is deterministic, then Nexp=1.
            Nbs = 100;                         % number of bootstrapped samples to compute standard error
            
            % collects 95% bootstrap ci's for estimands
            ci_sigma_u = zeros(2,n_u)*nan;
            ci_xi_u = zeros(2,n_u)*nan;                        
            ci_p_nea_u = zeros(2,n_u)*nan;
            
            % arrays for storing estimands
            param_save = zeros(2,n_u);
            pc = zeros(1,n_u)*nan;                                                            % collects estimated collision probability for each threshold
            ue_save = zeros(1,n_u)*nan;                                                       % collects estimated upper endpoint
            max_data = max(max(trans_data));                                                % largest observed value
            negL = @(par, exceed_data,u) -sum( log(gppdf(exceed_data,par(2),par(1),u)) );   % negative log likelihood fcn.
            logit = 1;                                                                      % set to plot logarithm of p_nea when magnitude of p_nea varies alot

            % limits for plots of empirical and model distribution functions
            xplot_lower = 0;
            xplot_upper = 1;
            n_eval_cdf = 1000;    

            for k=1:n_u
                
                param = est_par_w_shake(trans_data(:), init, U(k), shake_var);
                param_save(:,k) = param;
                p_u = sum(trans_data(:) > U(k))/length(data_matrix(:));
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

                    for j=1:Nbs

                        param_bs = bootstrap_est(trans_data, Nenc, init, U(k), shake_var);
                        p_u = length(exceed)/length(data);
                        sigma_sample(j) = param_bs(1); %#ok<*SAGROW>
                        xi_sample(j) = param_bs(2); 
                        p_nea_sample(j) = p_u*(1 - gpcdf(trans(0), param_bs(2), param_bs(1),U(k)) );

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
    plot(1:10   , param_save(2,:))
    title('xi_{est}')
    hold on
    if compute_ci == 1
        plot(1:10,ci_xi_u)
        thr_index = thr_autofind(ci_xi_u, param_save(2,:), n_u);
        plot(thr_index, param_save(2,thr_index),'*');
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

    if save_plot == 1
        savefig(sprintf('range_of_threshold_plots_DAFEA_exp_%d_id4',jj))
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

        % compute hit rate
        hit_rate = sum(pc_save_matrix>0,1)/500*100;

        % save results
        if select_trans == 1
            save(sprintf('hit_rate_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',data_type, select_trans, [], safety_level, Nenc),'hit_rate')
            save(sprintf('pc_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',      data_type, select_trans, [], safety_level, Nenc), 'pc_save_matrix')
            save(sprintf('p_c_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',     data_type, select_trans, [], safety_level, Nenc), 'p_c_save_matrix')
            save(sprintf('thr_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',     data_type, select_trans, [], safety_level, Nenc), 'thr_save_matrix')
            save(sprintf('param_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',   data_type, select_trans, [], safety_level, Nenc), 'param_save_matrix')
            save(sprintf('param_ci_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',data_type, select_trans, [], safety_level, Nenc), 'ci_xi_u_matrix')
            save(sprintf('p_exceed_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',data_type, select_trans, [], safety_level, Nenc), 'p_exceed_matrix')    

        elseif select_trans == 2
            save(sprintf('hit_rate_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',data_type, select_trans, p_ex*10, safety_level, Nenc),'hit_rate')
            save(sprintf('pc_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',      data_type, select_trans, p_ex*10, safety_level, Nenc), 'pc_save_matrix')
            save(sprintf('p_c_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',     data_type, select_trans, p_ex*10, safety_level, Nenc), 'p_c_save_matrix')
            save(sprintf('param_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',   data_type, select_trans, p_ex*10, safety_level, Nenc), 'param_save_matrix')
            save(sprintf('thr_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',     data_type, select_trans, p_ex*10, safety_level, Nenc), 'thr_save_matrix')
            save(sprintf('param_ci_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',data_type, select_trans, p_ex*10, safety_level, Nenc), 'ci_xi_u_matrix')
            save(sprintf('p_exceed_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',     data_type, select_trans, p_ex*10, safety_level, Nenc), 'p_exceed_matrix')    
        elseif select_trans == 3
            save(sprintf('hit_rate_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_samplesize_%d',data_type, select_trans, p_inv, d_inv, safety_level, Nenc),'hit_rate')
            save(sprintf('pc_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_samplesize_%d',      data_type, select_trans, p_inv, d_inv, safety_level, Nenc), 'pc_save_matrix')
            save(sprintf('p_c_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_samplesize_%d',     data_type, select_trans, p_inv, d_inv, safety_level, Nenc), 'p_c_save_matrix')
            save(sprintf('param__datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_samplesize_%d',  data_type, select_trans, p_inv, d_inv, safety_level, Nenc), 'param_save_matrix')
            save(sprintf('thr_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',        data_type, select_trans, p_inv, d_inv, safety_level, Nenc), 'thr_save_matrix')
            save(sprintf('param_ci_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',   data_type, select_trans, p_inv, d_inv, safety_level, Nenc), 'ci_xi_u_matrix')
            save(sprintf('p_exceed_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_samplesize_%d',   data_type, select_trans, p_inv, d_inv, safety_level, Nenc), 'p_exceed_matrix')    

        end

    end

end

