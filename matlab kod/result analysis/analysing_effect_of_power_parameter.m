%% analysing hit rate
hit_rate = get_data(1, 1, 2, 0.2, 1, 80, 6);
hit_rate = get_data(1, 1, 3, [4 3.5], 1, 80, 6);


plot(hit_rate)

%% analysing mean value
pc2 = get_data(3, 1, 2, 0.1, 1, 80, 6);
ci2 = compute_ci_meanest(pc2, 500);
pc4 = get_data(3, 1, 2, 0.2, 1, 80, 6);
ci4 = compute_ci_meanest(pc4, 500);
pc6 = get_data(3, 1, 2, 0.3, 1, 80, 6);
ci6 = compute_ci_meanest(pc6, 500);
pc8 = get_data(3, 1, 2, 0.4, 1, 80, 6);
ci8 = compute_ci_meanest(pc8, 500);
pc10 = get_data(3, 1, 2, 0.5, 1, 80, 6);
ci10 = compute_ci_meanest(pc10, 500);

load 'pc_nea_2mill_r_0_3.mat'
p_true = pc_nea;

clf
a = plot(mean(pc2), 'color' ,[0 .0 1]);
hold on
plot(ci2,'LineStyle', ':', 'color' ,[0 .0 1])
b = plot(mean(pc4), 'color' ,[0 .2 1]);
plot(ci4,'LineStyle', ':', 'color' ,[0 .2 1])
c = plot(mean(pc6),'color',[0 .4 1]);
plot(ci6,'LineStyle', ':', 'color' ,[0 .4 1])
d = plot(mean(pc8),'color',[0 .6 1]);
plot(ci8,'LineStyle', ':', 'color' ,[0 .6 1])
e = plot(mean(pc10), 'color' ,[0 .8 1]);
plot(ci10,'LineStyle', ':', 'color' ,[0 .8 1])

line([0,10], [p_true,p_true],'color','r')
title('mean value for different transformation parameters')
legend([a b c d e],'0.2','0.4','0.6','0.8','1.0')
ylim([0 3e-3])
xlabel('threshold number')
% we can see that there is clearly increase in bias as trans-par increases.
% specifically, the tendency to overestimate increases. In order to
% compensate for this increased bias, we have to use higher thresholds.
%% analysing Mean-Square-Error
pc2 = get_data(3, 1, 2, 0.2, 1);
pc4 = get_data(3, 1, 2, 0.4, 1);
pc6 = get_data(3, 1, 2, 0.6, 1);
pc8 = get_data(3, 1, 2, 0.8, 1);
pc10 = get_data(3, 1, 2, 1, 1);


rmse2 = sqrt(sum((pc2-p_true).^2)/500);
rmse4 = sqrt(sum((pc4-p_true).^2)/500);
rmse6 = sqrt(sum((pc6-p_true).^2)/500);
rmse8 = sqrt(sum((pc8-p_true).^2)/500);
rmse10 = sqrt(sum((pc10-p_true).^2)/500);
clf
plot(rmse2)
hold on
plot(rmse4)
plot(rmse6)
plot(rmse8)
plot(rmse10)
legend('0.2','0.4','0.6','0.8','1.0')
title('root-mean-square-error')
xlabel('threshold number')
%% comments on rmse
% As we might expect from the bias plots, the p=0.2 has the lowest and most
% consistent rmse. For the other parameters, we see a clear decline in bias
% as threshold increases.

%% analysing standard deviation
pc2 = get_data(3, 1, 2, 0.1, 1, 80, 6);
pc4 = get_data(3, 1, 2, 0.2, 1, 80, 6);
pc6 = get_data(3, 1, 2, 0.3, 1, 80, 6);
pc8 = get_data(3, 1, 2, 0.4, 1, 80, 6);
pc10 = get_data(3, 1, 2, 0.5, 1, 80, 6);


std2 = std(pc2);
std4 = std(pc4);
std6 = std(pc6);
std8 = std(pc8);
std10 = std(pc10);

clf
plot(std2)
hold on
plot(std4)
plot(std6)
plot(std8)
plot(std10)
legend('0.2','0.4','0.6','0.8','1.0')
title('standard deviation')
xlabel('threshold number')

%% comparing accuracy for different parameters

cut_off = .9;

pc2 = get_data(3, 1, 2, 0.2, 1);
pc4 = get_data(3, 1, 2, 0.4, 1);
pc6 = get_data(3, 1, 2, 0.6, 1);
pc8 = get_data(3, 1, 2, 0.8, 1);
pc10 = get_data(3, 1, 2, 1, 1);
pc_neg_ttc = get_data(3, 3, 1, [], 1);
pc_neg_stoch_ttc = get_data(3, 1, 1, [], 1);
pc_neg_min_dist = get_data(3, 6, 1, [], 1);

accuracy2 = accuracy_rate(pc2, p_true, cut_off);
accuracy4 = accuracy_rate(pc4, p_true, cut_off);
accuracy6 = accuracy_rate(pc6, p_true, cut_off);
accuracy8 = accuracy_rate(pc8, p_true, cut_off);
accuracy10 = accuracy_rate(pc10, p_true, cut_off);
accuracy_neg_ttc = accuracy_rate(pc_neg_ttc, p_true, cut_off);
accuracy_neg_stoch_ttc = accuracy_rate(pc_neg_stoch_ttc, p_true, cut_off);
accuracy_neg_min_dist = accuracy_rate(pc_neg_min_dist, p_true, cut_off);

ci2 = compute_ci_pest(accuracy2,500);
ci4 = compute_ci_pest(accuracy4,500);
ci6 = compute_ci_pest(accuracy6,500);
ci8 = compute_ci_pest(accuracy8,500);
ci10 = compute_ci_pest(accuracy10,500);
ci_neg_ttc = compute_ci_pest(accuracy_neg_ttc,500);
ci_neg_stoch_ttc = compute_ci_pest(accuracy_neg_stoch_ttc,500);
ci_neg_min_dist = compute_ci_pest(accuracy_neg_min_dist,500);


clf
a=plot(100*accuracy2,'color', [0 0 1]);
hold on
plot(100*ci2, ':','color', [0 0 1])
b=plot(100*accuracy4,'color', [0 .25 1]);
plot(100*ci4, ':','color','g')
c=plot(100*accuracy6,'color', [0 .5 1]);
plot(100*ci6, ':','color','g')
d=plot(100*accuracy8,'color', [0 .75 1]);
plot(100*ci8, ':','color','g')
e=plot(100*accuracy10,'color', [0 1 1]);
plot(100*ci10, ':','color','g')
f=plot(100*accuracy_neg_ttc,'color', 'black');
plot(100*ci_neg_ttc, ':','color','black')
g=plot(100*accuracy_neg_stoch_ttc,'color', 'red');
plot(100*ci_neg_stoch_ttc, ':','color','red')
h=plot(100*accuracy_neg_min_dist,'color', 'green');
plot(100*ci_neg_min_dist, ':','color','green')

legend([a b c d e f],'0.2','0.4','0.6','0.8','1.0', 'negated data')

%% comments on accuracy plots
% for p=0.2, the accuracy is fairly stable accross different thresholds,
% whereas for p=1.0 it starts out at zero and then peaks at the second to
% last threshold. This likely reflects the fact that for strong power
% parameter, the method strongly overestimates, and we need to use higher
% threshold to compensate for this bias. Note that for different values of
% pp that p=1.0 at its best performs better than p=0.2, and that the effect
% becomes more noticable as pp aproaches 1. It seems that the effect of
% increasing the transformation parameter (within the range studied) is 
% that the peak performance increases, but also the performance becomes
% more sensitive to threshold, with lower minimum performance. The plot
% suggests that the most appropriate power-parameter to use is 0.4, as it
% is insensitive to threshold while having approximately the same peak
% performance as the other values. In other worlds, it appears to strike
% the best balance between stability and peak performance.

% So, to answer the questions posed in the thesis: There does seem to be
% significant differences in performance between different transformation 
% parameters with regard to this measure of accuracy. For instance, p=0.4
% has significantly higher accuracy for pp>=5 than p=0.2, especially for 
% pp close to 1. in other words, if we want high probability of getting an
% estimate that is roughly correct in this sense, we would likely be better
% of using p=0.4 over p=0.2.

% The question now is whether or not we can use this information to see if
% we can come up with a method to find the optimal parameter before seeing
% the results, for instance by looking at stability plot or plots 
% comparing the distribution functions.

% There seems to be a tradeoff at work here. On one hand, we may use a
% transformation that has very little bias on one hand, but has a high
% chance of delivering zero/estimates. On the other hand, we can have a
% stronger transformation which is highly biased except for when the
% threshold is high, but which is much less likely to deliver
% zero/estimates. This seems to suggests the following hypothesis: If the stability
% plot is stable across all the thresholds, then we should select a higher
% threshold. If the stability plot is highly skewed, then we might do
% better by relaxing the transformation a little. In other words, we want
% to use a transformation that is stable, but not too stable. We want to
% formulate an algorithm that strikes a balance between stability and 

% note that the negated data has significantly lower success rate.
% Comparing the transformation of p_ex = 0.4 to negated transform, we see
% that there is a difference of approximately 30% in peak performance! At
% peak performance, we would get a good approximation about 4 times more
% often by using this transformation.

%% using automated threshold selection to select threshold
%% comparing accuracy, fail/success rates and bias
p_true = find_true_p(1);
success_matrix = zeros(1,5);
failure_matrix = zeros(1,5);
mean_matrix = zeros(1,5);
accuracy_matrix = zeros(1,5);
accuracy_ci_matrix = zeros(5,2);
cutoff = 0.5;


clf
par_range = [1.5 2 2.5 3 3.5 5];
for nn=1:5
    par = par_range(nn);

    thr_ind1 = zeros(1,500);
    pc_est1 = zeros(1,500);
    pc1 =      get_data(3, 1, 3, par, 1, 80, 6);
    xi1 = imag(get_data(4, 1, 3, par, 1, 80, 6));
    ci1 =      get_data(5, 1, 3, par, 1, 80, 6);

    thr_ind2 = zeros(1,500);
    pc_est2 = zeros(1,500);
    pc2 =      get_data(3, 1, 3, par, 2, 80, 6);
    xi2 = imag(get_data(4, 1, 3, par, 2, 80, 6));
    ci2 =      get_data(5, 1, 3, par, 2, 80, 6);

    for k=1:500
        thr_ind1(k) = thr_autofind(ci1{k}, xi1(k,:), 10);
        thr_ind2(k) = thr_autofind(ci2{k}, xi2(k,:), 10);
        pc_est1(k) = pc1(k,thr_ind1(k));
        pc_est2(k) = pc1(k,thr_ind2(k));
    end
    [succ_rate, fail_rate] = est_succ_fail_rate(pc_est1,pc_est2);
    
    
    accuracy_matrix(nn) = accuracy_rate(pc_est1, p_true, cutoff);
    accuracy_ci_matrix(nn,:) = compute_ci_pest(accuracy_matrix(nn),500);
    
    success_matrix(nn) = succ_rate;
    failure_matrix(nn) = fail_rate;
    mean_matrix(nn) = mean(pc_est1);
    
    clf
    plot(pc_est1)
    line([0 500], [mean(pc_est1) mean(pc_est1)])
    pause(.0)
end

exp_par = linspace(.2,1,5);
subplot(311)
plot(exp_par,success_matrix,'*')
hold on
plot(exp_par,failure_matrix,'*')
title('success and failure rates')
xlabel('transformation parameter')

subplot(312)
plot(exp_par, mean_matrix)
hold on
plot(exp_par,mean_matrix,'*')
plot([.2 1],[p_true p_true])
title('mean value of crash-frequency-estimates')
xlabel('transformation parameter')

subplot(313)
plot(exp_par, 100*accuracy_matrix)
hold on
plot(exp_par, 100*accuracy_ci_matrix,':','color', 'blue')
title(sprintf('accuracy plot with cutoff value 0.%d', cutoff*100))
xlabel('transformation parameter')
ylim([min(100*accuracy_matrix)*0.9, max(100*accuracy_matrix)*1.1])

%% comments on using automatic threshold selection
% We see that when using the automatic threshold selection, the successrate
% is quite low. 
%% PLots of success and failure rates
m =500;
plots = cell(5,2);
data_markers = data_marker();
data_type = 5;

par_range = [.1 .2 .3 .4 .5];
par_range = [1.5 3 5];% 1.5 2 2.5 3];
trans = 3;
npar = length(par_range);
leg_str = cell(1,npar); for i=1:npar; leg_str{i} = strcat('p=',num2str(par_range(i)) );end
leg_str{end+1}='no transform';
succ_fail_rate = cell(npar+1,2);
ci = cell(npar+1,2);

clf
for k=1:npar + 1
    if k == npar + 1
        trans = 1;
    else
        par = par_range(k); 
    end
    
    pc_safe1 = get_data(3, 5, trans, par, 1, 80,6);
    pc_safe2 = get_data(3, 5, trans, par, 2, 80,6);

    [succ_fail_rate{k,1}, succ_fail_rate{k,2}] = est_succ_fail_rate(pc_safe1, pc_safe2);

    ci{k,1} = compute_ci_pest(succ_fail_rate{k,1}, 500);
    ci{k,2} = compute_ci_pest(succ_fail_rate{k,2}, 500);
end

clf
subplot(211)
hold on
for k = 1:npar+1
    if k == npar+1
        col = 'r';
    else
        col = [0 (k-1)/k 1];
    end
    plots{k,1} = plot(succ_fail_rate{k,1}*100,data_markers{k}, 'color' ,col);
    plot(succ_fail_rate{k,1}*100, 'color' ,col);
    plot(ci{k,1}*100, ':', 'color' ,col)
    title('fail rate for different transformation parameters')
    xlabel('threshold index')
    legend([plots{1,1},plots{2,1},plots{3,1},plots{4,1},plots{5,1}], leg_str)
end

subplot(212)
hold on
for k = 1:npar+1
    if k == npar+1
        col = 'r';
    else
        col = [0 (k-1)/k 1];
    end
    plots{k,2} = plot(succ_fail_rate{k,2}*100, data_markers{k}, 'color' ,col);
    plot(succ_fail_rate{k,2}*100, 'color' ,col);
    plot(ci{k,2}*100, ':', 'color' ,col)
    xlabel('threshold index')
    legend([plots{1,2},plots{2,2},plots{3,2},plots{4,2},plots{5,2}], leg_str)
    title('fail rate for different transformation parameters')
end
%% comments on success/failure plots
% Increasing the power parameter seems to have the effect of increasing the
% success rate. This is likely due to the fact that, as we have observed,
% increasing the power parameter has the effect of decreasing the NZE-rate,
% which results in fewer indeterminate tests. From these plots we can
% conclude that if we want to maximize success rate and mimimize failure
% rate, we should use p_ex=1.0, and use the lowest threshold. It is likely
% that we are better of avoiding the retunr level to estimate.

%% Using return levels to determine which intersection is safer
m =10;
plots = cell(5,2);
data_markers = data_marker();
data_type = 5;

par_range = [.1 .2 .3 .4 .5];
par_range = [1.5 3 5];% 1.5 2 2.5 3];
trans = 3;
npar = length(par_range);
leg_str = cell(1,npar); for i=1:npar; leg_str{i} = strcat('p=',num2str(par_range(i)) );end
leg_str{end+1}='no transform';
succ_fail_rate = cell(npar+1,2);
ci = cell(npar+1,2);

clf
for k=1:npar + 1
    if k == npar + 1
        trans = 1;
    else
        p_ex = par_range(k); 
    end
    
    
    param1 = get_data(4, data_type, trans, p_ex, 1, 80, 6);
    xi_matrix1 = imag(param1);
    sigma_matrix1 = real(param1);
    u_matrix1 = get_data(6, data_type, trans, p_ex, 1, 80, 6);
    p_exceed1 = get_data(7, data_type, trans, p_ex, 1, 80, 6);
    x_m1 = return_level_m(m, sigma_matrix1, xi_matrix1, u_matrix1, p_exceed1);    
    
    param2 = get_data(4, data_type, trans, p_ex, 2, 80, 6);
    xi_matrix2 = imag(param2);
    sigma_matrix2 = real(param2);
    u_matrix2 = get_data(6, data_type, trans, p_ex, 2, 80, 6);
    p_exceed2 = get_data(7, data_type, trans, p_ex, 2, 80, 6);
    x_m2 = return_level_m(m, sigma_matrix2, xi_matrix2, u_matrix2, p_exceed2);

    [succ_fail_rate{k,1}, succ_fail_rate{k,2}] = est_succ_fail_rate(x_m1, x_m2);

    ci{k,1} = compute_ci_pest(succ_fail_rate{k,1}, 500);
    ci{k,2} = compute_ci_pest(succ_fail_rate{k,2}, 500);
end

clf
subplot(211)
hold on
for k = 1:npar+1
    if k == npar+1
        col = 'r';
    else
        col = [0 (k-1)/k 1];
    end
    plots{k,1} = plot(succ_fail_rate{k,1}*100,data_markers{k}, 'color' ,col);
    plot(succ_fail_rate{k,1}*100, 'color' ,col);
    plot(ci{k,1}*100, ':', 'color' ,col)
    title(sprintf('success rate for different transformation parameters, return period = %d',m))
    xlabel('threshold index')
    legend([plots{1,1},plots{2,1},plots{3,1},plots{4,1},plots{5,1}], leg_str)
end

subplot(212)
hold on
for k = 1:npar+1
    if k == npar+1
        col = 'r';
    else
        col = [0 (k-1)/k 1];
    end
    plots{k,2} = plot(succ_fail_rate{k,2}*100, data_markers{k}, 'color' ,col);
    plot(succ_fail_rate{k,2}*100, 'color' ,col);
    plot(ci{k,2}*100, ':', 'color' ,col)
    xlabel('threshold index')
    legend([plots{1,2},plots{2,2},plots{3,2},plots{4,2},plots{5,2}], leg_str)
    title(sprintf('fail rate for different transformation parameters, return period = %d',m))
end

%% comments on using return levels for infering safety levels
% well... apparantly its pretty fucking obvious which method should be used
% to make inferences about safety levels... Using p_exp=1, we manage to get
% 100% success rate for 5 out of 10 thresholds, using m=100.
%
% Three effects can be observed from the plots. First, we still see that the
% stronger transformations actually improves the ability to infer safety
% levels, as performance seems to grow with the strength of the
% transformation. Why this phenomena occurs is hard to say. It is possible
% that it is due to the fact that increasing the power of the
% transformation has the effect of lessening the effect of outliers,
% resulting in more stability, as the inference depends more heavily on
% observations that are more regular. Second, we see that the inference is
% more accurate when using lower thresholds. Third, we see that lowering the 
% return period m results in significantly improved estimates.
%
% These observations together suggets that performance improves when we
% are basing the inference on more data. In other words, we are letting the
% regular observations do most of the talking. This actually makes sense,
% since - given what we know about the simulation - there is no reason why
% basing our guess on non-estreme observations would not work well for
% safety order inference, since infering the safety order is equivalent to inferring
% which data set corresponds to the larger radius - an inference for which
% there is no need to focus on extreme events.
%
% It is important to note however that this experiment does not necesarily
% demonstrate that we will see these benefits with traffic data. What we are
% seing is that in order to infer which traffic location is safer for such
% simulated environments, we can look at the regular data, i.e. we are
% dealing with a type of data for which inference about the extremes based
% on the not-extreme data happens to work well. It is not clear if the same
% would be true for say real traffic data; it is not clear to which extent
% we can say something about the relative danger levels of two
% intersections when basing the inference mostly on regular non-extreme
% behaviour. Put in other words, it is not clear how extreme the events in
% traffic have to be before they start providing us with valuable
% information about collision probabilities. One cannot simply assume
% continuity accross the spectrum of encounter-severity the way we can for
% this simulated data.






