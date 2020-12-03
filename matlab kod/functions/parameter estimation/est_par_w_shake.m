function param = est_par_w_shake(data, init_val, u, shake_var)
% similar to est_par, but will check if algorithm is stuck on initial
% guess. If it is stuck, it will randomly modify the initial guess to
% attempt to get unstuck.

param =  est_par(data(:), init_val, u);

init_temp = init_val;
while_counter = 0;

while param == init_temp  
    
    while_counter = while_counter + 1;
    sigma0 = max( 0.1, init_val(1) + normrnd(0,shake_var^2) );
    xi0 = init_val(2) + normrnd(0,shake_var^2);
    init_temp = [sigma0, xi0];
    param =  est_par(data(:), init_temp, u);
    
    if while_counter == 30
        break
    end
end

end