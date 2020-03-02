%% analysing mean value
pc = get_data(3, 1, 1, [], 1);
ci2 = compute_ci_meanest(pc, 500);



p_true = 8.5333e-05;

clf
a = plot(mean(pc), 'color' ,[0 .0 1]);
hold on
plot(ci2,'LineStyle', ':', 'color' ,[0 .0 1])


line([0,10], [p_true,p_true],'color','r')
title('mean value for different transformation parameters')
xlabel('threshold number')

%% analysing Mean-Square-Error
pc2 = get_data(3, 1,1, [], 1);



rmse2 = sqrt(sum((pc2-p_true).^2)/500);

clf
plot(rmse2)
title('root-mean-square-error')
xlabel('threshold number')
%%
p_true = 8.5333e-05;
cut_off = 1;

pc_neg = get_data(3, 1, 1, [], 1);


accuracy_neg = accuracy_rate(pc_neg, p_true, cut_off);


ci_neg = compute_ci_pest(accuracy_neg,500);



clf
a=plot(100*accuracy2,'color', 'black');
hold on
plot(100*ci2, ':','color', 'black')




