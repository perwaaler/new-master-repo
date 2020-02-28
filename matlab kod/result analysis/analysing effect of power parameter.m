%% analysing hit rate
load hit_rate_datatype_1_trans_2_transpar_2_safetylevel_2_samplesize_334.mat
hit_rate2 = hit_rate
load hit_rate_datatype_1_trans_2_transpar_4_safetylevel_2_samplesize_334.mat
hit_rate4 = hit_rate
load hit_rate_datatype_1_trans_2_transpar_6_safetylevel_2_samplesize_334.mat
hit_rate6 = hit_rate
load hit_rate_datatype_1_trans_2_transpar_8_safetylevel_2_samplesize_334.mat
hit_rate8 = hit_rate
load hit_rate_datatype_1_trans_2_transpar_10_safetylevel_2_samplesize_334.mat
hit_rate10 = hit_rate

plot(hit_rate2)
hold on
plot(hit_rate4)
plot(hit_rate6)
plot(hit_rate8)
plot(hit_rate10)
%% analysing mean value
load p_c_datatype_1_trans_2_transpar_2_safetylevel_2_samplesize_334.mat
pc2 = p_c_save_matrix;
ci2 = compute_ci_meanest(p_c_save_matrix, 500)
load p_c_datatype_1_trans_2_transpar_4_safetylevel_2_samplesize_334.mat
pc4 = p_c_save_matrix;
ci4 = compute_ci_meanest(p_c_save_matrix, 500)
load p_c_datatype_1_trans_2_transpar_6_safetylevel_2_samplesize_334.mat
pc6 = p_c_save_matrix;
ci6 = compute_ci_meanest(p_c_save_matrix, 500)
load p_c_datatype_1_trans_2_transpar_8_safetylevel_2_samplesize_334.mat
pc8 = p_c_save_matrix;
ci8 = compute_ci_meanest(p_c_save_matrix, 500)
load p_c_datatype_1_trans_2_transpar_10_safetylevel_2_samplesize_334.mat
pc10 = p_c_save_matrix;
ci10 = compute_ci_meanest(p_c_save_matrix, 500)


p_true = 8.5333e-05;

clf
a = plot(mean(pc2), 'color' ,[0 .0 1])
hold on
plot(ci2,'LineStyle', ':', 'color' ,[0 .0 1])

b = plot(mean(pc4), 'color' ,[0 .2 1])
plot(ci4,'LineStyle', ':', 'color' ,[0 .2 1])

c = plot(mean(pc6),'color',[0 .4 1])
plot(ci6,'LineStyle', ':', 'color' ,[0 .4 1])

d = plot(mean(pc8),'color',[0 .6 1])
plot(ci8,'LineStyle', ':', 'color' ,[0 .6 1])

e = plot(mean(pc10), 'color' ,[0 .8 1])
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
load p_c_datatype_1_trans_2_transpar_2_safetylevel_2_samplesize_334.mat
pc2 = p_c_save_matrix;
load p_c_datatype_1_trans_2_transpar_4_safetylevel_2_samplesize_334.mat
pc4 = p_c_save_matrix;
load p_c_datatype_1_trans_2_transpar_6_safetylevel_2_samplesize_334.mat
pc6 = p_c_save_matrix;
load p_c_datatype_1_trans_2_transpar_8_safetylevel_2_samplesize_334.mat
pc8 = p_c_save_matrix;
load p_c_datatype_1_trans_2_transpar_10_safetylevel_2_samplesize_334.mat
pc10 = p_c_save_matrix;


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
load p_c_datatype_1_trans_2_transpar_2_safetylevel_2_samplesize_334.mat
pc2 = p_c_save_matrix;
load p_c_datatype_1_trans_2_transpar_4_safetylevel_2_samplesize_334.mat
pc4 = p_c_save_matrix;
load p_c_datatype_1_trans_2_transpar_6_safetylevel_2_samplesize_334.mat
pc6 = p_c_save_matrix;
load p_c_datatype_1_trans_2_transpar_8_safetylevel_2_samplesize_334.mat
pc8 = p_c_save_matrix;
load p_c_datatype_1_trans_2_transpar_10_safetylevel_2_samplesize_334.mat
pc10 = p_c_save_matrix;


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

%% analysing accuracy %%

pp = 0.9;

frac2 = pc2/p_true;
frac4 = pc4/p_true;
frac6 = pc6/p_true;
frac8 = pc8/p_true;
frac10 = pc10/p_true;


clf
plot(frac2(:,8))
hold on
%plot(frac4(:,1))
plot(frac6(:,8))
%plot(frac8(:,1))
%plot(frac10(:,1))
line([0,500], [1,1],'LineStyle','--','color','r')

accuracy2 = mean(abs(frac2 - 1)<pp);
accuracy4 = mean(abs(frac4 - 1)<pp);
accuracy6 = mean(abs(frac6 - 1)<pp)
accuracy8 = mean(abs(frac8 - 1)<pp)
accuracy10 = mean(abs(frac10 - 1)<pp)

ci2 = compute_ci_pest(accuracy2,500)
ci4 = compute_ci_pest(accuracy4,500)
ci6 = compute_ci_pest(accuracy6,500)
ci8 = compute_ci_pest(accuracy8,500)
ci10 = compute_ci_pest(accuracy10,500)


clf
a=plot(accuracy2,'r')
hold on
plot(ci2, ':','color','r')
b=plot(accuracy4,'b')
plot(ci4, ':','color','g')
c=plot(accuracy6,'m')
plot(ci6, ':','color','g')
d=plot(accuracy8,'color','black')
plot(ci8, ':','color','g')
e=plot(accuracy10,'g')
plot(ci10, ':','color','g')
legend([a b c d e],'0.2','0.4','0.6','0.8','1.0')

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
%% using automated threshold selection to select threshold
load p_c_datatype_1_trans_2_transpar_2_safetylevel_2_samplesize_334.mat
pc2 = p_c_save_matrix;
load param_datatype_1_trans_2_transpar_2_safetylevel_2_samplesize_334.mat
load param_ci_datatype_1_trans_2_transpar_2_safetylevel_2_samplesize_334.mat
xi_matrix = imag(param_save_matrix)
all_ci2 = ci_xi_u_matrix

thr_autofind(all_ci2{2,1},xi,n_thr)









