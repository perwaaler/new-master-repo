function out = get_data(variable, data_type, select_trans, p_ex, safety_level)
% function that loads the desired data-set. 
% set variable=1 for hit_rate
% set variable=2 for pc
% set variable=3 for p_c
% set variable=4 for param
% set variable=5 for param_ci
% if safety 


if variable == 1
    load(sprintf('hit_rate_datatype_%d_trans_%d_transpar_%d_safetylevel_%d',data_type, select_trans, p_ex*10, safety_level))
    out = hit_rate;
elseif variable == 2
    load(sprintf('pc_datatype_%d_trans_%d_transpar_%d_safetylevel_%d',data_type, select_trans, p_ex*10, safety_level)) %#ok<*LOAD>
    out = pc_save_matrix;
elseif variable == 3
    load(sprintf('p_c_datatype_%d_trans_%d_transpar_%d_safetylevel_%d',data_type, select_trans, p_ex*10, safety_level))
    out = p_c_save_matrix;
elseif variable == 4
    load(sprintf('param_datatype_%d_trans_%d_transpar_%d_safetylevel_%d',data_type, select_trans, p_ex*10, safety_level))
    out = param_save_matrix;
elseif variable == 5
    load(sprintf('param_ci_datatype_%d_trans_%d_transpar_%d_safetylevel_%d',data_type, select_trans, p_ex*10, safety_level))
    out = ci_xi_u_matrix;
elseif variable == 6
    load(sprintf('thr_datatype_%d_trans_%d_transpar_%d_safetylevel_%d',data_type, select_trans, p_ex*10, safety_level))
    out = thr_save_matrix;
elseif variable == 7
    load(sprintf('p_exceed_datatype_%d_trans_%d_transpar_%d_safetylevel_%d',data_type, select_trans, p_ex*10, safety_level))
    out = p_exceed_matrix;
end

end