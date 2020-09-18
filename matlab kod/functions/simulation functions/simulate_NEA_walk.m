function out = simulate_NEA_walk(A0,B0,stepsize, theta, driver_prop)
% simulates a random walk for points A and B, using initial conditions
% given by A0, B0, stepsize, theta, and driver_prop

save_walk = inf*ones(2,500);
counter = 1;
save_walk(:,1) = [A0,B0];

while imag(A0)<4 || real(B0)>-6
    newstep = take_NEA_step(A0,B0, stepsize, theta, driver_prop);
    A1 = newstep{1}(1);    
    B1 = newstep{1}(2); 
    stepsize = newstep{2};
    theta = newstep{3};
    counter = counter + 1;
    save_walk(:,counter) = [A1;B1];
    A0 = A1;
    B0 = B1;
end

out = save_walk;
end