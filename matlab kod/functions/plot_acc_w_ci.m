function a = plot_acc_w_ci(data_type, transform, trans_par, p_true, cut_off, color, data_mark)

pc = get_data(3, data_type, transform, trans_par, 1, 80, 6);
accuracy = accuracy_rate(pc, p_true, cut_off);
ci = compute_ci_pest(accuracy,500);

a=plot(100*accuracy, data_mark, 'color', color);
hold on
plot(100*accuracy,'color', color);
plot(100*ci,':','color', color)

end



