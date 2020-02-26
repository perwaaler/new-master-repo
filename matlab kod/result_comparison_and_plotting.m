%% compare hitrates of experiment 4 and 5
load hit_rate_ex_DAFEA_6_id4.mat
hit_rate_id4= hit_rate;
load hit_rate_ex_DAFEA_6_id5.mat
hit_rate_id5 = hit_rate;
clf
plot(hit_rate_id4);hold on;plot(hit_rate_id5)
legend('id4','id5')
title('percentage of non-zero-estimates useing DAFEA')
%% compare p_pure_coll (DAFEA and danger_FEA) from experiment 4 and 5, exp_par = 0.6
load p_c_ex_DAFEA_6_id4.mat
pc_pure_id4 = p_c_save_matrix;
load p_c_ex_DAFEA_6_id5.mat
pc_pure_id5 = p_c_save_matrix;
success_rate_DAFEA = sum(pc_pure_id5>pc_pure_id4)/500;
fail_rate_DAFEA = sum(pc_pure_id5<pc_pure_id4)/500

% compare p_pure_coll (danger_FEA) from experiment 4 and 5
load p_c_ex_danger_FEA_6_id4.mat
pc_pure_id4 = p_c_save_matrix;
load p_c_ex_danger_FEA_6_id5.mat
pc_pure_id5 = p_c_save_matrix;
success_rate_danger_FEA = sum(pc_pure_id5>pc_pure_id4)/500;
fail_rate_danger_FEA = sum(pc_pure_id5<pc_pure_id4)/500;

ci_success_DAFEA = success_rate_DAFEA' +  (sqrt( success_rate_DAFEA.*(1-success_rate_DAFEA)/500 )*1.96)'*[-1 1]
ci_success_danger_FEA = success_rate_danger_FEA' +  (sqrt( success_rate_danger_FEA.*(1-success_rate_danger_FEA)/500 )*1.96)'*[-1 1]
ci_fail_DAFEA = fail_rate_DAFEA' +  (sqrt( fail_rate_DAFEA.*(1-fail_rate_DAFEA)/500 )*1.96)'*[-1 1]
ci_fail_danger_FEA = fail_rate_danger_FEA' +  (sqrt( fail_rate_danger_FEA.*(1-fail_rate_danger_FEA)/500 )*1.96)'*[-1 1]


clf
subplot(211)
a = plot(success_rate_DAFEA,'r');hold on
plot(ci_success_DAFEA,'--','color','r')
b = plot(success_rate_danger_FEA,'b')
plot(ci_success_danger_FEA,'--','color','b')
legend([a b],'DAFEA','dangerFEA')
title('success rate, p_exp = 0.6')

subplot(212)
a = plot(fail_rate_DAFEA,'r');hold on
plot(ci_fail_DAFEA,'--','color','r')
b = plot(fail_rate_danger_FEA,'b')
plot(ci_fail_danger_FEA,'--','color','b')
legend([a b],'DAFEA','dangerFEA')
title('fail rate, p_exp = 0.6')
%% comparing p_exp=0.6 against p_exp=0.7 in ability to guess most dangerous intersection

load p_c_ex_DAFEA_6_id4.mat
pc_pure_6_id4 = p_c_save_matrix;
load p_c_ex_DAFEA_6_id5.mat
pc_pure_6_id5 = p_c_save_matrix;
fails_DAFEA_6 = pc_pure_6_id5<pc_pure_6_id4;
fail_rate_6 = sum(fails_DAFEA_6)/500;

load p_c_ex_DAFEA_7_id4.mat
pc_pure_7_id4 = p_c_save_matrix;
load p_c_ex_DAFEA_7_id5.mat
pc_pure_7_id5 = p_c_save_matrix;
fails_DAFEA_7 = pc_pure_7_id5<pc_pure_7_id4;
fail_rate_7 = sum(fails_DAFEA_7)/500;
%%% make confidence intervals %%%
ci_6 = fail_rate_6' +  (sqrt( fail_rate_6.*(1-fail_rate_6)/500 )*1.96)'*[-1 1]
ci_7 = fail_rate_7' + (sqrt( fail_rate_7.*(1-fail_rate_7)/500 )*1.96)'*[-1 1]

clf; plot(fail_rate_6,'r');hold on;plot(ci_6,'--','color','r') 
plot(fail_rate_7,'b');hold on;plot(ci_7,'--','color','b')
legend('pexp=0.6','pexp=0.7')
title('failure rate for pexp=0.6 and pexp=0.7')

%% comparing X against DAFEA in ability to guess most dangerous intersection

% DAFEA
load p_c_ex_DAFEA_6_id4.mat
pc_pure_DAFEA_id4 = p_c_save_matrix;
load p_c_ex_DAFEA_6_id5.mat
pc_pure_DAFEA_id5 = p_c_save_matrix;
fail_rate_DAFEA = sum(pc_pure_DAFEA_id5<pc_pure_DAFEA_id4)/500;
success_rate_DAFEA = sum(pc_pure_DAFEA_id5>pc_pure_DAFEA_id4)/500;

% X
load p_c_ex_X_6_id4.mat
pc_pure_X_id4 = p_c_save_matrix;
load p_c_ex_X_6_id5.mat
pc_pure_X_id5 = p_c_save_matrix;
fail_rate_X = sum(pc_pure_X_id5<pc_pure_X_id4)/500;
success_rate_X = sum(pc_pure_X_id5>pc_pure_X_id4)/500;

%%% make confidence intervals %%%
ci_fail_DAFEA = fail_rate_DAFEA' + sqrt( fail_rate_DAFEA.*(1-fail_rate_DAFEA)/500 )'*[-1 1]*1.96
ci_fail_X = fail_rate_X' + sqrt( fail_rate_X.*(1-fail_rate_X)/500 )'*[-1 1]*1.96
ci_success_DAFEA = success_rate_DAFEA' + sqrt( success_rate_DAFEA.*(1-success_rate_DAFEA)/500 )'*[-1 1]*1.96
ci_success_X = success_rate_X' + sqrt( success_rate_X.*(1-success_rate_X)/500 )'*[-1 1]*1.96

clf; 
subplot(211)
a = plot(fail_rate_DAFEA,'r');hold on
plot(ci_fail_DAFEA,'--','color','r')
b = plot(fail_rate_X,'b')
plot(ci_fail_X,'--','color','b')
legend([a,b],'DAFEA','X')
title('failure rate for X and DAFEA, using T(x)=exp(-0.6x)')

subplot(212)
a = plot(success_rate_DAFEA,'r');hold on;
plot(ci_success_DAFEA,'--','color','r');
b = plot(success_rate_X,'b');
plot(ci_success_X,'--','color','b');
legend([a,b],'DAFEA','X')
title('success rate for X and DAFEA, using T(x)=exp(-0.6x)')

%% comparing failure rates when using neg DAFEA and X with upper endpoints for comparison

load ue_neg_DAFEA_id4.mat
ue_id4 = ue_save_matrix;
load ue_neg_DAFEA_id5.mat
ue_id5 = ue_save_matrix;
fail_rate_ue_DAFEA = sum(ue_id5<ue_id4)/500;
success_rate_ue_DAFEA = sum(ue_id5>ue_id4)/500;

load ue_neg_X_id4.mat
ue_id4 = ue_save_matrix;
load ue_neg_X_id5.mat
ue_id5 = ue_save_matrix;
fail_rate_ue_X = sum(ue_id5<ue_id4)/500;
success_rate_ue_X = sum(ue_id5>ue_id4)/500;


ci_fail_DAFEA =  fail_rate_ue_DAFEA' + (sqrt( fail_rate_ue_DAFEA.*(1-fail_rate_ue_DAFEA)/500 )*1.96)'*[-1 1]
ci_success_DAFEA =  success_rate_ue_DAFEA' + (sqrt( success_rate_ue_DAFEA.*(1-success_rate_ue_DAFEA)/500 )*1.96)'*[-1 1]
ci_fail_X =  fail_rate_ue_X' + (sqrt( fail_rate_ue_X.*(1-fail_rate_ue_X)/500 )*1.96)'*[-1 1]
ci_success_X =  success_rate_ue_X' + (sqrt( success_rate_ue_X.*(1-success_rate_ue_X)/500 )*1.96)'*[-1 1]

clf; 
subplot(211)
plot(fail_rate_ue_DAFEA,'r');hold on; plot(ci_fail_DAFEA,'--','color','r')
plot(fail_rate_ue_X,'b');hold on; plot(ci_fail_X,'--','color','b')
title('failure rate using upper-endpoint for comparison')
subplot(212)
plot(success_rate_ue_DAFEA,'r');hold on; plot(ci_success_DAFEA,'--','color','r')
plot(success_rate_ue_X,'b');hold on; plot(ci_success_X,'--','color','b')
title('success rate using upper-endpoint for comparison')
%% comparison in success/fail rate between transformed and untransformed data


% DAFEA transformed
load p_c_ex_DAFEA_6_id4.mat
pc_pure_DAFEA_id4 = p_c_save_matrix;
load p_c_ex_DAFEA_6_id5.mat
pc_pure_DAFEA_id5 = p_c_save_matrix;
fail_rate_DAFEA = sum(pc_pure_DAFEA_id5<pc_pure_DAFEA_id4)/500;
success_rate_DAFEA = sum(pc_pure_DAFEA_id5>pc_pure_DAFEA_id4)/500;

% DAFEA negated
load p_c_neg_DAFEA_id4.mat
pc_pure_DAFEA_neg_id4 = p_c_save_matrix;
load p_c_neg_DAFEA_id5.mat
pc_pure_DAFEA_neg_id5 = p_c_save_matrix;
fail_rate_DAFEA_neg = sum(pc_pure_DAFEA_neg_id5<pc_pure_DAFEA_neg_id4)/500;
success_rate_DAFEA_neg = sum(pc_pure_DAFEA_neg_id5>pc_pure_DAFEA_neg_id4)/500;

%%% make confidence intervals %%%
ci_fail_DAFEA = fail_rate_DAFEA' + sqrt( fail_rate_DAFEA.*(1-fail_rate_DAFEA)/500 )'*[-1 1]*1.96
ci_fail_DAFEA_neg = fail_rate_DAFEA_neg' + sqrt( fail_rate_DAFEA_neg.*(1-fail_rate_DAFEA_neg)/500 )'*[-1 1]*1.96
ci_success_DAFEA = success_rate_DAFEA' + sqrt( success_rate_DAFEA.*(1-success_rate_DAFEA)/500 )'*[-1 1]*1.96
ci_success_DAFEA_neg = success_rate_DAFEA_neg' + sqrt( success_rate_DAFEA_neg.*(1-success_rate_DAFEA_neg)/500 )'*[-1 1]*1.96

clf; 
subplot(211)
a = plot(fail_rate_DAFEA,'r');hold on
plot(ci_fail_DAFEA,'--','color','r')
b = plot(fail_rate_DAFEA_neg,'b')
plot(ci_fail_DAFEA_neg,'--','color','b')
legend([a,b],'DAFEA transformed','DAFEA negated')
title('fail rate for DAFEA transformed and negated, using T(x)=exp(-0.6x)')

subplot(212)
a = plot(success_rate_DAFEA,'r');hold on;
plot(ci_success_DAFEA,'--','color','r');
b = plot(success_rate_DAFEA_neg,'b');
plot(ci_success_DAFEA_neg,'--','color','b');
legend([a,b],'DAFEA transformed','DAFEA negated')
title('success rate for DAFEA transformed and negated, using T(x)=exp(-0.6x)')
%% comparing accuracy
spill = 0.9
p_emp = 8.5333e-05
load p_c_ex_DAFEA_6_id4.mat
pc_pure_DAFEA_id4 = p_c_save_matrix;
fail_rate_DAFEA = abs((1-pc_pure_DAFEA_id4/p_emp))<spill
accuracy_DAFEA_exp6 = sum(fail_rate_DAFEA)/500
ci_fail_DAFEA_ex6 = accuracy_DAFEA_exp6' + sqrt( accuracy_DAFEA_exp6.*(1-accuracy_DAFEA_exp6)/500 )'*[-1 1]*1.96

load p_c_ex_DAFEA_6_id4.mat
pc_pure_DAFEA_id4 = p_c_save_matrix;
fail_rate_DAFEA = abs((1-pc_pure_DAFEA_id4/8.8000e-05))<spill
accuracy_DAFEA_exp7 = sum(fail_rate_DAFEA)/500
load p_c_ex_X_6_id4.mat
ci_fail_DAFEA_ex7 = accuracy_DAFEA_exp7' + sqrt( accuracy_DAFEA_exp7.*(1-accuracy_DAFEA_exp7)/500 )'*[-1 1]*1.96

load p_c_ex_X_6_id4.mat
pc_pure_X_id4 = p_c_save_matrix;
fail_rate_X = abs((1-pc_pure_X_id4/p_emp))<spill
accuracy_X_exp6 = sum(fail_rate_X)/500
ci_fail_X_ex6 = accuracy_X_exp6' + sqrt( accuracy_X_exp6.*(1-accuracy_X_exp6)/500 )'*[-1 1]*1.96


load p_c_ex_danger_FEA_6_id4.mat
pc_pure_danger_FEA_id4 = p_c_save_matrix;
success_danger_FEA = abs((1-pc_pure_danger_FEA_id4/p_emp))<spill
accuracy_danger_FEA_exp6 = sum(success_danger_FEA)/500
ci_fail_danger_FEA_ex6 = accuracy_danger_FEA_exp6' + sqrt( accuracy_danger_FEA_exp6.*(1-accuracy_danger_FEA_exp6)/500 )'*[-1 1]*1.96


clf;
a = plot(accuracy_DAFEA_exp7,'r'); hold on
plot(ci_fail_DAFEA_ex7,'--','color','r')
%plot(accuracy_X_exp6,'b');
b = plot(accuracy_danger_FEA_exp6,'b');
%plot(ci_fail_X_ex6,'--','color','b')
plot(ci_fail_danger_FEA_ex6,'--','color','b')
title('probability to deviate less than 80% from true value')
legend([a b], 'stoch. TTC, transformed', 'sep. dist., transformed')
xlabel('threshold number')
%% mean value plots
load p_c_ex_DAFEA_6_id4.mat
pc_pure_DAFEA_id4 = p_c_save_matrix;
mean_val_DAFEA = mean(pc_pure_DAFEA_id4)
sdev_DAFEA = std(pc_pure_DAFEA_id4)/sqrt(500)

load p_c_ex_X_6_id4.mat
pc_pure_X_id4 = p_c_save_matrix;
mean_val_X = mean(pc_pure_X_id4)
sdev_X = std(pc_pure_X_id4)/sqrt(500)


load p_c_ex_danger_FEA_6_id4.mat
pc_pure_danger_FEA_id4 = p_c_save_matrix;
mean_val_danger_FEA = mean(pc_pure_danger_FEA_id4)
sdev_danger_FEA = std(pc_pure_danger_FEA_id4)/sqrt(500)

ci_mean_val_DAFEA = mean_val_DAFEA' + sdev_DAFEA'*[-1 1]*1.96
ci_mean_val_X = mean_val_X' + sdev_X'*[-1 1]*1.96
ci_mean_val_danger_FEA = mean_val_danger_FEA' + sdev_danger_FEA'*[-1 1]*1.96
clf
a = plot(mean_val_DAFEA,'r'); hold on; 
plot(ci_mean_val_DAFEA,'--','color','r')
b = plot(mean_val_X,'b'); hold on; 
plot(ci_mean_val_X,'--','color','b')
c = plot(mean_val_danger_FEA,'m'); hold on; 
plot(ci_mean_val_danger_FEA,'--','color','m')
legend([a b c], 'DAFEA, T(x) = exp(-0.6x)','X','dangerFEA')

% we can see in these plots that the price we pay for the imrovement in
% reduction of non-zero estimates comes at the cost of over-estimation.
% Danger FEA has significantly lower and more correct expectation value,
% but it also estimates the collision probability to be zero quite often.
% Hence being correct on average does not mean that it is a good estimator.

