function [row, col] = find_sim_sit(pos0, speed0, theta0, historyA, tolerance)

% tolerance = [0.15, 0.02, 0.2];
% pos0 = -1 - 1i;
% speed0 = 0.16;

tol_pos   = tolerance(1);
tol_speed = tolerance(2);
tol_theta = tolerance(3);

save_pos_A   = historyA{1};
save_step_A  = historyA{2};
save_theta_A = historyA{3};

NN = size(save_step_A,1);
MM = size(save_step_A,2);

xx = 1:NN;
XX = sqrt((save_pos_A - pos0).*conj(save_pos_A - pos0));
[minval,time] = min(XX,[],2);
enc_ind = xx(minval < tol_pos); % trajectories that pass through the disk centered at AA0 with radius tol_pos
time = time(minval<tol_pos);          % times at which the circle was passed through
coordinates = sub2ind([NN MM], enc_ind, time');
XX(sub2ind(size(save_step_A), enc_ind, time'));

% time and enc_ind gives the indeces for the encounters and time during
% those encounters that the condition was satisfied

speeds = save_step_A(coordinates);
passes_speed_test = find(abs(speeds-speed0)<tol_speed);
coordinates = coordinates(passes_speed_test);

thetas = save_theta_A(coordinates);
passes_theta_test = find(abs(thetas-theta0)<tol_theta);

coordinates = coordinates(passes_theta_test);
[row,col] = ind2sub(size(save_step_A),coordinates);
% %%
% iter=5
% walk = save_pos_A(row(iter),col(iter):end)
% walk = walk(save_pos_A(row(iter),col(iter):end)<inf)
% plot(walk)
% %%
end


%% check that code does the right thing!
% sqrt((save_pos_A(coordinates) - pos0).*conj(save_pos_A(coordinates) - pos0))
% abs(save_step_A(coordinates) - speed0)
% abs(save_theta_A(coordinates) - theta0)