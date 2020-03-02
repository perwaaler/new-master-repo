% 1=DAFEA 2=X 3=ttc_FEA 4=ttc_min ttc_min_EA=5 danger_FEA=6 dist_min_EA=7 dist_min=8

load 'pc_nea_2mill_r_0_3.mat'
load 'pc_ea_2mill_r_0_3.mat'

cut_off = 0.9;
pc_tot = pc_nea + pc_ea;
p_exp = 0.9;
pc_min_ttc = get_data(3, 3, 2, 0.3, 1, 80, 6);

pc_split = get_data(3, 5,2, 0.2, 1, 80, 6) + get_data(3, 3, 2, 0.1, 1, 80, 6);

p_combined = 1/2*(pc_split + pc_min_ttc)

accuracy = accuracy_rate(pc_min_ttc, pc_tot, cut_off);
accuracy_split = accuracy_rate(pc_split, pc_tot, cut_off);
accuracy_combined = accuracy_rate(p_combined, pc_tot, cut_off);

ci_min_ttc = compute_ci_pest(accuracy,500);
ci_split = compute_ci_pest(accuracy_split,500);
ci_combined = compute_ci_pest(accuracy_combined,500);

clf
a=plot(100*accuracy,'o','color', [0 0 1]);
hold on
plot(100*accuracy,'color', [0 0 1]);
plot(100*ci_min_ttc,':','color',[0 0 1])

b=plot(100*accuracy_split,'o','color', [0 1 1]);
plot(100*accuracy_split,'color', [0 1 1]);
plot(100*ci_split,':','color',[0 1 1])

c = plot(100*accuracy_combined,'o','color', 'red');
plot(100*accuracy_combined,'color', 'red');
plot(100*ci_combined,':','color','red')


legend([a b], 'minimum ttc', 'data splitting method, ttc')
%% using only 

% 1=DAFEA 2=X 3=ttc_FEA 4=ttc_min ttc_min_EA=5 danger_FEA=6 dist_min_EA=7 dist_min=8

load 'pc_nea_2mill_r_0_3.mat'
load 'pc_ea_2mill_r_0_3.mat'

cut_off = 0.95;
pc_tot = pc_nea + pc_ea;
p_exp = 0.9;
pc_min_ttc = get_data(3, 3, 1, [], 1, 80, 6);

pc_split = get_data(3, 1,1, [], 1, 80, 6) + get_data(3, 5, 1, [], 1, 80, 6);

p_combined = 1/2*(pc_split + pc_min_ttc)

accuracy = accuracy_rate(pc_min_ttc, pc_tot, cut_off);
accuracy_split = accuracy_rate(pc_split, pc_tot, cut_off);
accuracy_combined = accuracy_rate(p_combined, pc_tot, cut_off);

ci_min_ttc = compute_ci_pest(accuracy,500);
ci_split = compute_ci_pest(accuracy_split,500);
ci_combined = compute_ci_pest(accuracy_combined,500);

clf
a=plot(100*accuracy,'o','color', [0 0 1]);
hold on
plot(100*accuracy,'color', [0 0 1]);
plot(100*ci_min_ttc,':','color',[0 0 1])

b=plot(100*accuracy_split,'o','color', [0 1 1]);
plot(100*accuracy_split,'color', [0 1 1]);
plot(100*ci_split,':','color',[0 1 1])

c = plot(100*accuracy_combined,'o','color', 'red');
plot(100*accuracy_combined,'color', 'red');
plot(100*ci_combined,':','color','red')


legend([a b], 'minimum ttc', 'data splitting method, ttc')







