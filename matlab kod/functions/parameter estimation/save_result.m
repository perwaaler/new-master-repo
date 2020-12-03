function save_result(variables, data_type, select_trans, par, safety_level, up_frac, lo_frac)
% saves variables
hit_rate =          variables{1};
pc_save_matrix =    variables{2}; %#ok<*NASGU>
p_c_save_matrix =   variables{3};
thr_save_matrix =   variables{4};
param_save_matrix = variables{5};
ci_xi_u_matrix =    variables{6};
p_exceed_matrix =   variables{7};

if select_trans == 1
    save(sprintf('hit_rate_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, [], safety_level, up_frac*100, lo_frac*100), 'hit_rate')
    save(sprintf(      'pc_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, [], safety_level, up_frac*100, lo_frac*100), 'pc_save_matrix') %#ok<*USENS>
    save(sprintf(     'p_c_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, [], safety_level, up_frac*100, lo_frac*100), 'p_c_save_matrix')
    save(sprintf(     'thr_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, [], safety_level, up_frac*100, lo_frac*100), 'thr_save_matrix')
    save(sprintf(   'param_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, [], safety_level, up_frac*100, lo_frac*100), 'param_save_matrix')
    save(sprintf('param_ci_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, [], safety_level, up_frac*100, lo_frac*100), 'ci_xi_u_matrix')
    save(sprintf('p_exceed_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, [], safety_level, up_frac*100, lo_frac*100), 'p_exceed_matrix')    

elseif select_trans == 2
    save(sprintf('hit_rate_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, par*100, safety_level, up_frac*100, lo_frac*100), 'hit_rate')
    save(sprintf(      'pc_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, par*100, safety_level, up_frac*100, lo_frac*100), 'pc_save_matrix')
    save(sprintf(     'p_c_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, par*100, safety_level, up_frac*100, lo_frac*100), 'p_c_save_matrix')
    save(sprintf(   'param_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, par*100, safety_level, up_frac*100, lo_frac*100), 'param_save_matrix')
    save(sprintf(     'thr_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, par*100, safety_level, up_frac*100, lo_frac*100), 'thr_save_matrix')
    save(sprintf('param_ci_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, par*100, safety_level, up_frac*100, lo_frac*100), 'ci_xi_u_matrix')
    save(sprintf('p_exceed_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',    data_type, select_trans, par*100, safety_level, up_frac*100, lo_frac*100), 'p_exceed_matrix')    

elseif select_trans == 3
    save(sprintf('hit_rate_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d', data_type, select_trans, par(1)*10, par(2)*10, safety_level, up_frac*100, lo_frac*100), 'hit_rate')
    save(sprintf('pc_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d',       data_type, select_trans, par(1)*10, par(2)*10, safety_level, up_frac*100, lo_frac*100), 'pc_save_matrix')
    save(sprintf('p_c_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d',      data_type, select_trans, par(1)*10, par(2)*10, safety_level, up_frac*100, lo_frac*100), 'p_c_save_matrix')
    save(sprintf('param_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d',   data_type, select_trans, par(1)*10, par(2)*10, safety_level, up_frac*100, lo_frac*100), 'param_save_matrix')
    save(sprintf('thr_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d',      data_type, select_trans, par(1)*10, par(2)*10, safety_level, up_frac*100, lo_frac*100), 'thr_save_matrix')
    save(sprintf('param_ci_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d', data_type, select_trans, par(1)*10, par(2)*10, safety_level, up_frac*100, lo_frac*100), 'ci_xi_u_matrix')
    save(sprintf('p_exceed_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d', data_type, select_trans, par(1)*10, par(2)*10, safety_level, up_frac*100, lo_frac*100), 'p_exceed_matrix') 
end

end
    