function p_true = find_true_p(data_type)
if data_type == 1 || data_type == 3 || data_type == 6
    load 'pc_nea_2mill_r_0_3.mat'
    p_true = pc_nea;
elseif data_type == 5 || data_type == 7 
    load 'pc_ea_2mill_r_0_3.mat'
    p_true = pc_ea;
else
    load 'pc_2mill_r_0_3.mat'
    p_true = pc_tot;
end