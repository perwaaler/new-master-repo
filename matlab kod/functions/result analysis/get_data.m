function out = get_data(variable, data_type, select_trans, trans_par, safety_level, up_frac, lo_frac)
% function that loads the desired data-set. 
% set variable=1 for hit_rate
% set variable=2 for pc
% set variable=3 for p_c
% set variable=4 for param
% set variable=5 for param_ci
% if safety 

if length(trans_par)==1 && select_trans==3
    trans_par=[trans_par 3.5];
end
if select_trans == 1
    trans_par=[];
end
    
if select_trans==1 || select_trans==2
if variable == 1
    load(sprintf('hit_rate_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par*100, safety_level, up_frac, lo_frac))
    out = hit_rate;
elseif variable == 2
    load(      sprintf('pc_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par*100, safety_level, up_frac, lo_frac)) %#ok<*LOAD>
    out = pc_save_matrix;
elseif variable == 3
    load(     sprintf('p_c_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par*100, safety_level, up_frac, lo_frac))
    out = p_c_save_matrix;
elseif variable == 4
    load(   sprintf('param_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par*100, safety_level, up_frac, lo_frac))
    out = param_save_matrix;
elseif variable == 5
    load(sprintf('param_ci_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par*100, safety_level, up_frac, lo_frac))
    out = ci_xi_u_matrix;
elseif variable == 6
    load(     sprintf('thr_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par*100, safety_level, up_frac, lo_frac))
    out = thr_save_matrix;
elseif variable == 7
    load(sprintf('p_exceed_datatype_%d_trans_%d_transpar_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par*100, safety_level, up_frac, lo_frac))
    out = p_exceed_matrix;
end
end

if select_trans==3
if variable == 1
    load(sprintf('hit_rate_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par(1)*10, trans_par(2)*10, safety_level, up_frac, lo_frac))
    out = hit_rate;
elseif variable == 2
    load(      sprintf('pc_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par(1)*10, trans_par(2)*10, safety_level, up_frac, lo_frac))
    out = pc_save_matrix;
elseif variable == 3
    load(     sprintf('p_c_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par(1)*10, trans_par(2)*10, safety_level, up_frac, lo_frac))
    out = p_c_save_matrix;
elseif variable == 4
    load(   sprintf('param_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par(1)*10, trans_par(2)*10, safety_level, up_frac, lo_frac))
    out = param_save_matrix;
elseif variable == 5
    load(sprintf('param_ci_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par(1)*10, trans_par(2)*10, safety_level, up_frac, lo_frac))
    out = ci_xi_u_matrix;
elseif variable == 6
    load(     sprintf('thr_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par(1)*10, trans_par(2)*10, safety_level, up_frac, lo_frac))
    out = thr_save_matrix;
elseif variable == 7
    load(sprintf('p_exceed_datatype_%d_trans_%d_transpar_%d_%d_safetylevel_%d_u_frac_%d_l_frac_%d',data_type, select_trans, trans_par(1)*10, trans_par(2)*10, safety_level, up_frac, lo_frac))
    out = p_exceed_matrix;
end
end



end