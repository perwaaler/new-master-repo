function trans_data = transform(data_matrix, select_trans, trans_par)
% function that transforms the data using the desired monotone transform.

if select_trans == 1
    trans_data = -data_matrix;
    
elseif select_trans == 2
    trans_data = exp( -trans_par*(data_matrix - 2));
    
elseif select_trans == 3
    trans_data = 1./(data_matrix + trans_par(1)).^ trans_par(2);
end
end