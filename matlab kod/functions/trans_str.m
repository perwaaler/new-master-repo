function str = trans_str(trans_ind)
% function that translates transformation index into a string

trans_string = cell(3,1);
trans_string{1} = 'no trans.';
trans_string{2} = 'exp.';
trans_string{3} = 'inv.';
str = trans_string{trans_ind};

end