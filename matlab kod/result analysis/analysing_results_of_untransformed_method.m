
%% plot estimated probabilities
clf
pc1 = get_data(3, 1, 3, 1, 1, 80, 6);
pc2 = get_data(3, 1, 3, 4.5, 1, 80, 6);
thr = 9;
p_true=find_true_p(1);
subplot(211)
plot(pc1(pc1(:,thr)>0,thr),'.')
hold on
plot(ones(1,sum(pc1(:,thr)>0))*p_true*0.5,'g')
plot(ones(1,sum(pc1(:,thr)>0))*p_true)
plot(ones(1,sum(pc1(:,thr)>0))*p_true*1.5,'g')
% mean value conditional on NZE
plot(ones(1,sum(pc1(:,thr)>0))*mean(pc1(pc1(:,thr)>0,thr)),'black');

subplot(212)
plot(pc2(pc2(:,thr)>0,thr),'.')
hold on
plot(ones(1,sum(pc2(:,thr)>0))*p_true*0.5,'g')
plot(ones(1,sum(pc2(:,thr)>0))*p_true)
plot(ones(1,sum(pc2(:,thr)>0))*p_true*1.5,'g')
plot(ones(1,sum(pc2(:,thr)>0))*mean(pc2(pc2(:,thr)>0,thr)),'black');


%% plot probabilities 
clf
pc_neg_stoch_ttc_exp2 = get_data(3, 1, 2, 0.25, 1, 80, 6);
pc_neg_stoch_ttc_inv30 = get_data(3, 1, 3, [4.5 3.5], 1, 80, 6);

p_true=find_true_p(1);

subplot(211)
plot(pc_neg_stoch_ttc_inv30(:,6)/p_true,'.')
hold on
plot(ones(1,500)*0.5,'r')
plot(ones(1,500)*1.5,'r')
subplot(212)
plot(pc_neg_stoch_ttc_exp2(:,6)/p_true,'.')
hold on
plot(ones(1,500)*0.5,'r')
plot(ones(1,500)*1.5,'r')
%% analysing mean value for untransformed methods
pc_stochttc =      get_data(3, 1, 1, [], 1, 80, 6);
pc_ttc     =       get_data(3, 2, 1, [], 1, 80, 6);
pc_stoch_ttc_exp = get_data(3, 1, 2, 0.1, 1, 80, 6);

ci_stoch_ttc = compute_ci_meanest(pc_stochttc, 500);
ci_ttc = compute_ci_meanest(pc_ttc, 500);
ci_stoch_ttc_exp = compute_ci_meanest(pc_stoch_ttc_exp, 500);



load 'pc_nea_2mill_r_0_3.mat'
p_true = pc_nea;

clf
a = plot(mean(pc_stochttc), 'color' ,[0 .0 1]);
hold on
plot(ci_stoch_ttc,'LineStyle', ':', 'color' ,[0 .0 1])

b = plot(mean(pc_ttc), 'color' ,'red');
plot(ci_ttc,'LineStyle', ':', 'color' ,'red')

a = plot(mean(pc_stoch_ttc_exp), 'color' ,[0 .0 1]);
hold on
plot(ci_stoch_ttc_exp,'LineStyle', ':', 'color' ,[0 .0 1])

line([0,10], [p_true,p_true],'color','r')
title('mean value for different transformation parameters')
xlabel('threshold number')

%% analysing standard deviation
pc1 = get_data(3, 1, 2, 0.1, 1, 80, 6);
pc2 = get_data(3, 1, 2, 0.2, 1, 80, 6);
pc3 = get_data(3, 1, 2, 0.3, 1, 80, 6);
pc4 = get_data(3, 1, 2, 0.4, 1, 80, 6);
pc5 = get_data(3, 1, 2, 0.5, 1, 80, 6);

pcinv10 = get_data(3, 1, 3, 1, 1, 80, 6);
pcinv15 = get_data(3, 1, 3, 1.5 , 1, 80, 6);
pcinv20 = get_data(3, 1, 3,  2 , 1, 80, 6);
pcinv25 = get_data(3, 1, 3, 2.5, 1, 80, 6);
pcinv30 = get_data(3, 1, 3, 3 , 1, 80, 6);
pcinv35 = get_data(3, 1, 3, 3.5, 1, 80, 6);
pcinv40 = get_data(3, 1, 3, 4 , 1, 80, 6);
pcinv45 = get_data(3, 1, 3, 4.5, 1, 80, 6);


std1 = std(pc1);
std2 = std(pc2);
std3 = std(pc3);
std4 = std(pc4);
std5 = std(pc5);
stdinv10 = std(pcinv10);
stdinv15 = std(pcinv15);
stdinv20 = std(pcinv20);
stdinv25 = std(pcinv25);
stdinv30 = std(pcinv30);
stdinv35 = std(pcinv35);
stdinv40 = std(pcinv40);
stdinv45 = std(pcinv45);

clf
subplot(211)

plot(std1, 'color', [0 0 1])
hold on
plot(std2, 'color', [0 .2 1])
plot(std3, 'color', [0 .4 1])
plot(std4, 'color', [0 .6 1])
plot(std5, 'color', [0 .8 1])

legend('p = 0.1','p = 0.2','p = 0.3','p = 0.4','p = 0.5')
title('standard deviation, exp. transform')
xlabel('threshold number')
subplot(212)
plot(stdinv10, 'color', [0 0 1])
hold on
plot(stdinv15, 'color', [0 .1 1])
plot(stdinv20, 'color', [0 .2 1])
plot(stdinv25, 'color', [0 .3 1])
plot(stdinv30, 'color', [0 .4 1])
plot(stdinv35, 'color', [0 .5 1])
plot(stdinv40, 'color', [0 .6 1])
plot(stdinv45, 'color', [0 .7 1])
legend('p = 1','p = 1.5','p = 2.0','p = 2.5','p = 3','p = 3.5','p = 4.0','p = 4.5')
title('standard deviation, inv. transform')
xlabel('threshold number')

%% analysing Mean-Square-Error
pc_stochttc = get_data(3, 3, 1, [3 3.5], 1, 80, 6);
%load 'pc_ea_2mill_r_0_3.mat'
p_true = pc_nea;

rmse2 = sqrt(sum((pc_stochttc-p_true).^2)/500);

hold on
% clf
plot(rmse2)
title('root-mean-square-error')
xlabel('threshold number')
%% comparison of severity measures at FEA
% load 'pc_nea_2mill_r_0_3.mat'
load 'pc_nea_1mill.mat'
p_true = pc_nea;
data_markers = data_marker();
sev_measure = cell(3,1);
sev_measure{1} = 'stoch. TTC';
sev_measure{2} = 'TTC';
sev_measure{3} = 'min. dist.';

trans_string = cell(3,1);
trans_string{1} = [];
trans_string{2} = 'exponential';
trans_string{3} = 'inverse';

cut_off = .9;

clf
a = plot_acc_w_ci(1 , 1, [3], p_true, cut_off, 'r', data_markers{1});
b = plot_acc_w_ci(3, 1, [.1], p_true, cut_off, 'b', data_markers{2});
c = plot_acc_w_ci(6, 1, [0.1], p_true, cut_off, 'm', data_markers{3});
title(sprintf('accuracy plots for severity measures, cut-off %s',num2str(cut_off)))
legend([a b c], sev_measure{1},sev_measure{2},sev_measure{3})

%  why does min dist perform so well?

%% comparing accuracy plots for exponential transform
% 1=DAFEA 2=X 3=ttc_FEA 4=ttc_min ttc_min_EA=5 danger_FEA=6 dist_min_EA=7 dist_min=8

cut_off = .5;
data_type = 1;
tran = 2;
par_range = [.1 .2 .3 .4 .5];
% par_range = [.1 .15 .2 .25 .3 .4 .5];  % ttc parameters
% par_range = [.02 .06 .1 .2 .3 .4 ]; % mindist parameters
par_range = [0.2];
% par_range = [1 1.5 2 2.5 3 3.5 4 4.5];

sev_ind = sev_measure(data_type);
%p_true = find_true_p(data_type);
n_ex = length(par_range);
data_markers = data_marker();


pc_stochttc_exp =        cell(1,n_ex);
accuracy_stoch_ttc_exp = cell(1,n_ex);
ci_stoch_ttc_exp =       cell(1,n_ex);
legend_labels =          cell(n_ex, 1);
color_spec =             linspace(0,1,n_ex);
plot_vec =               zeros(1, n_ex);

clf
hold on
for i=1:n_ex
    
    pc_stochttc_exp{i} = get_data(3, data_type, tran, par_range(i), 1, 80, 6);
    accuracy_stoch_ttc_exp{i} = accuracy_rate(pc_stochttc_exp{i}, p_true, cut_off);
    ci_stoch_ttc_exp{i} = compute_ci_pest(accuracy_stoch_ttc_exp{i},500);


    plot_vec(i)=plot(100*accuracy_stoch_ttc_exp{i}, data_markers{i},'color', [0 color_spec(i) 1]);
    plot(100*accuracy_stoch_ttc_exp{i},'color', [0 color_spec(i) 1]);
    plot(100*ci_stoch_ttc_exp{i},':','color',[0 color_spec(i) 1])
    
    if i==n_ex
        legend_labels{i,1} = 'no trans.';
    else
        legend_labels{i,1} = sprintf('p = %s', num2str(par_range(i)));
    end
end

a = plot_acc_w_ci(data_type, 1, [], p_true, cut_off, 'r', data_markers{8}, .8, 0.06);
%plot_acc_w_ci(data_type, transform, trans_par, p_true, cut_off, color, data_mark)
%b = plot_acc_w_ci(3, 1, [], p_true, cut_off, 'r', data_markers{9});
%c = plot_acc_w_ci(1, 3, [4.5 3.5], p_true, cut_off, 'r', data_markers{9});







legend_labels{end + 1} = sevme_str(sev_ind);
% labels{9} = 'ttc';

title(sprintf('accuracy plots for %s, %s transform, cut-off = 0.%d',sevme_str(sev_ind), trans_str(tran), cut_off*10))
xlabel('threshold index')
% ylim([0,17])
legend([plot_vec a], legend_labels)%labels)
%% accuracy inv transform
% 1=DAFEA 2=X 3=ttc_FEA 4=ttc_min ttc_min_EA=5 danger_FEA=6 dist_min_EA=7 dist_min=8



cut_off = .5;
data_type = 2;
tran = 3;
sev_ind = sev_measure(data_type);

p_true = find_true_p(data_type);

par_range = [ 1.5 2 3 3.5 4 4.5];
par_range = [.1 .3 .7 .9 1.2 1.4]; % par range for min dist
par_range = [1.5 3 4.5 5 5.5];
% par_range = [1 1.5 2 2.5 3 3.5 4 4.5];
n_ex = length(par_range);
data_markers = data_marker();


pc_stochttc_exp = cell(1,n_ex);
accuracy_stoch_ttc_exp = cell(1,n_ex);
ci_stoch_ttc_exp = cell(1,n_ex);
legend_labels = cell(n_ex, 1);
color_spec = linspace(0,1,n_ex);
plot_vec = zeros(1, n_ex);

clf
hold on
for i=1:n_ex

    pc_stochttc_exp{i} = get_data(3, data_type, tran, [par_range(i),3.5], 1, 80, 6);
    accuracy_stoch_ttc_exp{i} = accuracy_rate(pc_stochttc_exp{i}, p_true, cut_off);
    ci_stoch_ttc_exp{i} = compute_ci_pest(accuracy_stoch_ttc_exp{i},500);


    plot_vec(i)=plot(100*accuracy_stoch_ttc_exp{i}, data_markers{i},'color', [0 color_spec(i) 1]);
    plot(100*accuracy_stoch_ttc_exp{i},'color', [0 color_spec(i) 1]);
    plot(100*ci_stoch_ttc_exp{i},':','color',[0 color_spec(i) 1])

    legend_labels{i,1} = sprintf('p = %s', num2str(par_range(i)));
end

%a = plot_acc_w_ci(data_type, 3, [.1 3.5], p_true, cut_off, 'r', data_markers{end-1},.85, .06);
a = plot_acc_w_ci(data_type, 1, [], p_true, cut_off, 'r', data_markers{end-1}, .80, .06);

legend_labels{n_ex+1} = 'no trans.';

title(sprintf('accuracy plots for %s, %s transform, cut-off = %s', sevme_str(sev_ind),trans_str(tran) ,num2str(cut_off)) );
%ylim([0 25])
legend([plot_vec a], legend_labels)

%% comparison of peak performance between exponential and inverse transform
load 'pc_nea_2mill_r_0_3.mat'
p_true = pc_nea;

cut_off = .5;
clf
%plot_acc_w_ci(data_type, transform, trans_par, p_true, cut_off, color, data_mark)
hold on
a = plot_acc_w_ci(1, 2, 0.25, p_true, cut_off, 'r', 's');
b = plot_acc_w_ci(1, 3, [4 3.5], p_true, cut_off, 'b', 's');
legend([a b], 'stoch. TTC, inv, p = 4.0', 'stoch. TTC, exp, p = 0.25)')
title(sprintf('comparison of peak performance, cut-off = %s',num2str(cut_off)))

%% analysing mean value of inverse transform


cut_off = .5;
data_type = 1;
tran = 3;

sev_ind = sev_measure(data_type);
p_true = find_true_p(data_type);
par_range = [1 1.5 2 2.5 3 3.5 4 4.5]; % par range for stoch ttc
% par_range = [.1 .3 .7 .9 1.2 1.4];     % mindist parameters, inv
% par_range = [.02 .06 .1 .2 .3 .4 ];     % mindist parameters, exp
n_ex = length(par_range);
data_markers = data_marker();

pc_stochttc_exp = cell(1,n_ex);
ci_stoch_ttc_exp = cell(1,n_ex);
legend_labels = cell(n_ex, 1);
color_spec = linspace(0,1,n_ex);
plot_vec = zeros(1, n_ex);

clf
hold on
for i=1:n_ex

    pc_stochttc_exp{i} = get_data(3, data_type, tran, par_range(i), 1, 80, 6);
    accuracy_stoch_ttc_exp{i} = mean(pc_stochttc_exp{i});
    ci_stoch_ttc_exp{i} = compute_ci_meanest(pc_stochttc_exp{i},500);


    plot_vec(i)=plot(accuracy_stoch_ttc_exp{i}, data_markers{i},'color', [0 color_spec(i) 1]);
    plot(accuracy_stoch_ttc_exp{i},'color', [0 color_spec(i) 1]);
    plot(ci_stoch_ttc_exp{i},':','color',[0 color_spec(i) 1])

    legend_labels{i,1} = sprintf('p = %s',num2str(par_range(i)));
end

pc_stochttc = get_data(3, data_type, 1, [], 1, 80, 6);
accuracy_stoch_ttc = mean(pc_stochttc);
ci_stoch_ttc = compute_ci_meanest(pc_stochttc,500);

a = plot(accuracy_stoch_ttc, data_markers{i+2},'color', 'r');
plot(accuracy_stoch_ttc,'color', 'r');
plot(ci_stoch_ttc,':','color','r')

b = plot(ones(1,10)*p_true,'black');

legend_labels{n_ex+1} = 'no trans.';
legend_labels{n_ex+2} = pc_str(coll_type(data_type));

title(sprintf('Expected value for %s transform, cut-off = %s',trans_str(tran),num2str(cut_off)))
xlabel('threshold index')
legend([plot_vec a b], legend_labels)

%% analysing mean value of exp transform

p_true = find_true_p(data_type);

cut_off = .5;
data_type = 1;
tran = 2;
sev_ind = sev_measure(data_type);

par_range = [.1 .15 .2 .25 .3 .4 .5];     % stoch ttc parameters
% par_range = [.02 .06 .1 .2 .3 .4 ];     % mindist parameters
n_ex = length(par_range);
data_markers = data_marker();

pc_stochttc_exp = cell(1,n_ex);
ci_stoch_ttc_exp = cell(1,n_ex);
legend_labels = cell(n_ex, 1);
color_spec = linspace(0,1,n_ex);
plot_vec = zeros(1, n_ex);

clf
hold on
for i=1:n_ex

    pc_stochttc_exp{i} = get_data(3, data_type, tran, par_range(i), 1, 80, 6);
    accuracy_stoch_ttc_exp{i} = mean(pc_stochttc_exp{i});
    ci_stoch_ttc_exp{i} = compute_ci_meanest(pc_stochttc_exp{i},500);


    plot_vec(i)=plot(accuracy_stoch_ttc_exp{i}, data_markers{i},'color', [0 color_spec(i) 1]);
    plot(accuracy_stoch_ttc_exp{i},'color', [0 color_spec(i) 1]);
    plot(ci_stoch_ttc_exp{i},':','color',[0 color_spec(i) 1])

    legend_labels{i,1} = sprintf('p = %s',num2str(par_range(i)));
end

pc_stochttc = get_data(3, data_type, 1, [], 1, 80, 6);
accuracy_stoch_ttc = mean(pc_stochttc);
ci_stoch_ttc = compute_ci_meanest(pc_stochttc,500);

a = plot(accuracy_stoch_ttc, data_markers{i+2},'color', 'r');
plot(accuracy_stoch_ttc,'color', 'r');
plot(ci_stoch_ttc,':','color','r')

ylim([0 1e-3])
b = plot(ones(1,10)*p_true,'black');

legend_labels{n_ex+1} = 'no transform';
legend_labels{n_ex+2} = 'P(NEA,C)';

title(sprintf('Expected value for exp. transform, cut-off %s',num2str(cut_off)))
xlabel('threshold index')
legend([plot_vec a b], legend_labels)
% we can see that there is clearly increase in bias as trans-par increases.
% specifically, the tendency to overestimate increases. In order to
% compensate for this increased bias, we have to use higher thresholds.
%% analysing hitrate
p_true = find_true_p(data_type);

cut_off = .5;
data_type = 1;
tran = 2;
sev_ind = sev_measure(data_type);

par_range = [.1 .15 .2 .25 .3 .4 .5];     % stoch ttc parameters
% par_range = [.02 .06 .1 .2 .3 .4 ];     % mindist parameters
% par_range = [1 1.5 2.0 2.5 3 3.5 4 4.5];
n_ex = length(par_range);
data_markers = data_marker();

hit_rate = cell(1,n_ex);
hit_rate_ci = cell(1,n_ex);
legend_labels = cell(n_ex, 1);
color_spec = linspace(0,1,n_ex);
plot_vec = zeros(1, n_ex);

clf
hold on
for i=1:n_ex

    hit_rate{i} = get_data(1, data_type, tran, par_range(i), 1, 80, 6);
    hit_rate_ci{i} = compute_ci_pest(hit_rate{i},500);


    plot_vec(i)=plot(hit_rate{i}, data_markers{i},'color', [0 color_spec(i) 1]);
    plot(hit_rate{i},'color', [0 color_spec(i) 1]);
    %plot(hit_rate_ci{i},':','color',[0 color_spec(i) 1])

    legend_labels{i,1} = sprintf('p = %s' ,num2str(par_range(i)));
end

hit_rate_notrans = get_data(1, data_type, 1, [], 1, 80, 6);
% ci_stoch_ttc = compute_ci_pest(hit_rate_notrans,500);

a = plot(hit_rate_notrans, data_markers{i+2},'color', 'r');
plot(hit_rate_notrans,'color', 'r');
% plot(ci_stoch_ttc,':','color','r')


legend_labels{n_ex+1} = 'no trans.';


title(sprintf('P(NZE) for %s transform, cut-off %s',trans_str(tran) ,num2str(cut_off)))
xlabel('threshold index')
legend([plot_vec a], legend_labels)
%% analysing hitrate
hitratettcstoch = get_data(1, 1, 1, [], 1, 80, 6)
hitratettc = get_data(1, 3, 1, [], 1, 80, 6)
hitratemindist = get_data(1, 6, 1, [], 1, 85, 6)

clf
hold on
plot(hitratettc)
plot(hitratettcstoch)
plot(hitratemindist)
%% analysing conditional hitrate
p_true=find_true_p(2);
pc = get_data(3, 1, 2, 0.1, 1, 80, 6);

n_good_approx = accuracy_rate(pc, p_true, cut_off)*500
n_NZE = sum(pc>0)

p_good_approx_cond = n_good_approx./n_NZE
%%







