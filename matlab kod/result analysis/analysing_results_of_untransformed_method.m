%% plot probabilities 
clf
pc_neg_stoch_ttc_exp2 = get_data(3, 1, 2, 0.2, 1, 80, 6);
pc_neg_stoch_ttc_inv30 = get_data(3, 1, 3, [3 3.5], 1, 80, 6);

load 'pc_ea_2mill_r_0_3.mat'
p_true = pc_nea;
plot(pc_neg_stoch_ttc_inv30(:,1)/p_true,'.')

hold on
plot(ones(1,500)*0.1)
plot(ones(1,500)*1.9)
%% analysing mean value
pc_stochttc = get_data(3, 1, 1, [], 1, 80, 6);
pc_ttc =       get_data(3, 3, 1, [], 1, 80, 6);
pc_stoch_ttc_exp = get_data(3, 1, 2, 0.1, 1, 80, 6);

ci_stoch_ttc = compute_ci_meanest(pc_stochttc, 500);
ci_ttc = compute_ci_meanest(pc_ttc, 500)
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

%% analysing Mean-Square-Error
pc_stochttc = get_data(3, 1, 3, [3 3.5], 1, 80, 6);
load 'pc_ea_2mill_r_0_3.mat'
p_true = pc_nea;

rmse2 = sqrt(sum((pc_stochttc-p_true).^2)/500);

clf
plot(rmse2)
title('root-mean-square-error')
xlabel('threshold number')
%% comparing accuracy plots
load 'pc_nea_2mill_r_0_3.mat'
p_true = pc_nea;

cut_off = .5;

par_range = [.1 .15 .2 .25 .3 .4 .5];
n_ex = length(par_range);
data_markers = cell(1,n_ex);
data_markers{1} = 'o';
data_markers{2} = '*';
data_markers{3} = '+';
data_markers{4} = 'x';
data_markers{5} = 'p';
data_markers{6} = 'd';
data_markers{7} = 's';
data_markers{8} = 'h';
data_markers{9} = '^';
data_markers{9} = 'v';


pc_stochttc_exp = cell(1,n_ex);
accuracy_stoch_ttc_exp = cell(1,n_ex);
ci_stoch_ttc_exp = cell(1,n_ex);
labels = cell(n_ex, 1);
color_spec = linspace(0,1,n_ex);
plot_vec = zeros(1, n_ex);

clf
hold on
for i=1:n_ex

    pc_stochttc_exp{i} = get_data(3, 1, 2, par_range(i), 1, 80, 6);
    accuracy_stoch_ttc_exp{i} = accuracy_rate(pc_stochttc_exp{i}, p_true, cut_off);
    ci_stoch_ttc_exp{i} = compute_ci_pest(accuracy_stoch_ttc_exp{i},500);


    plot_vec(i)=plot(100*accuracy_stoch_ttc_exp{i}, data_markers{i},'color', [0 color_spec(i) 1]);
    plot(100*accuracy_stoch_ttc_exp{i},'color', [0 color_spec(i) 1]);
    plot(100*ci_stoch_ttc_exp{i},':','color',[0 color_spec(i) 1])

    labels{i,1} = sprintf('stoch ttc, p = 0.%d',par_range(i)*100);
end

a = plot_acc_w_ci(1, 1, [], p_true, cut_off, 'r', data_markers{8});
%plot_acc_w_ci(data_type, transform, trans_par, p_true, cut_off, color, data_mark)
%b = plot_acc_w_ci(3, 1, [], p_true, cut_off, 'r', data_markers{9});
%c = plot_acc_w_ci(1, 3, [4.5 3.5], p_true, cut_off, 'r', data_markers{9});


labels{8} = 'stoch ttc';
labels{9} = 'ttc';

title(sprintf('accuracy plots for exponential transform, cutoff = 0.%d',cut_off*10))

legend([plot_vec a], labels)
%% accuracy inv transform
%
load 'pc_nea_2mill_r_0_3.mat'
p_true = pc_nea;

cut_off = .5;

par_range = [1 1.5 2 2.5 3 3.5 4 4.5];
n_ex = length(par_range);
data_markers = cell(1,n_ex);
data_markers{1} = 'o';
data_markers{2} = '*';
data_markers{3} = '+';
data_markers{4} = 'x';
data_markers{5} = 'p';
data_markers{6} = 'd';
data_markers{7} = 's';
data_markers{8} = 'h';
data_markers{9} = '^';
data_markers{9} = 'v';

pc_stochttc_exp = cell(1,n_ex);
accuracy_stoch_ttc_exp = cell(1,n_ex);
ci_stoch_ttc_exp = cell(1,n_ex);
labels = cell(n_ex, 1);
color_spec = linspace(0,1,n_ex);
plot_vec = zeros(1, n_ex);

clf
hold on
for i=1:n_ex

    pc_stochttc_exp{i} = get_data(3, 1, 3, [par_range(i),3.5], 1, 80, 6);
    accuracy_stoch_ttc_exp{i} = accuracy_rate(pc_stochttc_exp{i}, p_true, cut_off);
    ci_stoch_ttc_exp{i} = compute_ci_pest(accuracy_stoch_ttc_exp{i},500);


    plot_vec(i)=plot(100*accuracy_stoch_ttc_exp{i}, data_markers{i},'color', [0 color_spec(i) 1]);
    plot(100*accuracy_stoch_ttc_exp{i},'color', [0 color_spec(i) 1]);
    plot(100*ci_stoch_ttc_exp{i},':','color',[0 color_spec(i) 1])

    labels{i,1} = sprintf('stoch ttc, %d.%d',floor(par_range(i)), decpart(par_range(i)));;
end

a = plot_acc_w_ci(1, 1, [], p_true, cut_off, 'r', data_markers{end-1});



labels{n_ex+1} = 'stoch ttc';


title(sprintf('accuracy plots for inverse transform, cutoff = 0.%d',cut_off*10))

legend([plot_vec a], labels)

%% comparison of peak performance between exponential and inverse transform
load 'pc_nea_2mill_r_0_3.mat'
p_true = pc_nea;

cut_off = .5;
clf
%plot_acc_w_ci(data_type, transform, trans_par, p_true, cut_off, color, data_mark)
hold on
a = plot_acc_w_ci(1, 2, 0.25, p_true, cut_off, 'r', 's');
b = plot_acc_w_ci(1, 3, [4 3.5], p_true, cut_off, 'b', 's');
legend([a b], 'stoch. ttc, inv. trans., p = 4.0', 'stoch. ttc, exp, p = 0.25)')
title('comparison between best performing transformations')

%% analysing mean value of inverse transform

load 'pc_nea_2mill_r_0_3.mat'
p_true = pc_nea;

par_range = [1 1.5 2 2.5 3 3.5 4 4.5];
n_ex = length(par_range);
data_markers = cell(1,n_ex);
data_markers{1} = 'o';
data_markers{2} = '*';
data_markers{3} = '+';
data_markers{4} = 'x';
data_markers{5} = 'p';
data_markers{6} = 'd';
data_markers{7} = 's';
data_markers{8} = 'h';
data_markers{9} = '^';
data_markers{9} = 'v';

pc_stochttc_exp = cell(1,n_ex);
ci_stoch_ttc_exp = cell(1,n_ex);
labels = cell(n_ex, 1);
color_spec = linspace(0,1,n_ex);
plot_vec = zeros(1, n_ex);

clf
hold on
for i=1:n_ex

    pc_stochttc_exp{i} = get_data(3, 1, 3, [par_range(i),3.5], 1, 80, 6);
    accuracy_stoch_ttc_exp{i} = mean(pc_stochttc_exp{i});
    ci_stoch_ttc_exp{i} = compute_ci_meanest(pc_stochttc_exp{i},500);


    plot_vec(i)=plot(accuracy_stoch_ttc_exp{i}, data_markers{i},'color', [0 color_spec(i) 1]);
    plot(accuracy_stoch_ttc_exp{i},'color', [0 color_spec(i) 1]);
    plot(ci_stoch_ttc_exp{i},':','color',[0 color_spec(i) 1])

    labels{i,1} = sprintf('stoch ttc, %d.%d',floor(par_range(i)), decpart(par_range(i)));
end
a = plot(ones(1,10)*p_true,'r');

labels{n_ex+1} = 'pc_{nea} true';

title('mean value for inverse transform')

legend([plot_vec a], labels)

%% analysing mean value of exp transform




load 'pc_nea_2mill_r_0_3.mat'
p_true = pc_nea;

par_range = [.1 .15 .2 .25 .3 .4 .5];
n_ex = length(par_range);
data_markers = cell(1,n_ex);
data_markers{1} = 'o';
data_markers{2} = '*';
data_markers{3} = '+';
data_markers{4} = 'x';
data_markers{5} = 'p';
data_markers{6} = 'd';
data_markers{7} = 's';
data_markers{8} = 'h';
data_markers{9} = '^';
data_markers{9} = 'v';

pc_stochttc_exp = cell(1,n_ex);
ci_stoch_ttc_exp = cell(1,n_ex);
labels = cell(n_ex, 1);
color_spec = linspace(0,1,n_ex);
plot_vec = zeros(1, n_ex);

clf
hold on
for i=1:n_ex

    pc_stochttc_exp{i} = get_data(3, 1, 2, par_range(i), 1, 80, 6);
    accuracy_stoch_ttc_exp{i} = mean(pc_stochttc_exp{i});
    ci_stoch_ttc_exp{i} = compute_ci_meanest(pc_stochttc_exp{i},500);


    plot_vec(i)=plot(accuracy_stoch_ttc_exp{i}, data_markers{i},'color', [0 color_spec(i) 1]);
    plot(accuracy_stoch_ttc_exp{i},'color', [0 color_spec(i) 1]);
    plot(ci_stoch_ttc_exp{i},':','color',[0 color_spec(i) 1])

    labels{i,1} = sprintf('stoch ttc, p = %d.%d',floor(par_range(i)), decpart(par_range(i)));
end
ylim([0 2e-3])
a = plot(ones(1,10)*p_true,'r');

labels{n_ex+1} = 'pc_{nea} true';

title('mean value for exponential transform')

legend([plot_vec a], labels)
% we can see that there is clearly increase in bias as trans-par increases.
% specifically, the tendency to overestimate increases. In order to
% compensate for this increased bias, we have to use higher thresholds.


