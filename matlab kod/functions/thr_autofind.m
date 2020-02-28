function thr_opt = thr_autofind(conf_int,xi,n_thr)
% Takes the the matrix containing confidence
% intervals as columns, a vector xi containing the estimates and the
% number of thresholds n_est for which estimates are computed. Returns the
% index of the threshold with the lowest cost fcn value.

thr_opt = zeros(1,n_thr);
for i = 1:n_thr-1
    zum = 0;
    for j = i+1:10
        %zum = zum + triangle_fcn(xi(i) - xi(j), abs(conf_int(1,j)-conf_int(2,j))/2 )
        zum = zum + (abs(xi(i) - xi(j))/abs(conf_int(1,j)-conf_int(2,j))/2)^1.1;
        %zum = zum + triangle_fcn(xi(i) - xi(j), abs(conf_int(1,j)-conf_int(2,j))/2, 2);
    end
    thr_opt(i) = (mean(zum)+0.2)*(1 + 1.2*(i/n_thr)^2);
end
thr_opt %#ok<NOPRT>
[~,I] = min(thr_opt(1:end-1));
thr_opt = I;
end