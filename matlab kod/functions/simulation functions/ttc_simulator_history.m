function out = ttc_simulator_history(A0,B0, stepsize, theta, r, driver_prop, save_pos_A, row, col, iter, plot_setting)
% out = ttc_simulator(A0,B0,stepsize0,min_stepsize,sigma_step, r, thetavar)
% This function generates a random walk starting at position A0 and
% B0 respectively, and then computes Time To Collision. TTC is set to 1000
% if collision never occurs.
% A0 and B0 given current positions, and theta describe current direction.
% driver_prop = {stepsize0, var_step, var_theta, stability_fac, min_stepsize, theta_mod_par, speed_mod_par};



stepcounter = 0;
ttc = 1i;
D = norm(A0-B0);

% if they have already collided, run:
if D < 2*r
    % get pre-collision positions based on arguments
    B1 = B0 + stepsize(2)*exp(1i*theta(2));
    A1 = A0 - stepsize(1)*exp(1i*theta(1));
    Adiff = A1-A0;
    Bdiff = B1-B0;
    eta = real(conj(Adiff - Bdiff)*(A0-B0))/norm(Adiff-Bdiff)^2;
    mu = (norm(A0 - B0)^2 - 4*r^2)/norm(Adiff - Bdiff)^2;
    ttc = -eta + [-1 1]*sqrt(eta^2-mu);
    ttc = -max(ttc); % ttc will now be negative, since time must be reversed to find time of collision
    out = ttc;
else
    walk_id = row(iter);
    walk = save_pos_A(walk_id,col(iter):end);
    walk = walk(save_pos_A(row(iter),col(iter):end)<inf);
    while real(A0) < real(B0) + 0.9 && imag(A0)<imag(B0)+0.9
                
            % generate new step
            newstep = take_NEA_step(A0, B0, stepsize, theta, driver_prop);
            A1 = walk(stepcounter+2);
            B1 = newstep{2};
            stepsize = newstep{3};
            theta = newstep{4};

            if norm(A1-B1)<2*r
                Adiff = A1-A0;
                Bdiff = B1-B0;
                eta = real(conj(Adiff - Bdiff)*(A0-B0))/norm(Adiff-Bdiff)^2;
                mu = (norm(A0 - B0)^2 - 4*r^2)/norm(Adiff - Bdiff)^2;
                ttc = -eta + [-1 1]*sqrt(eta^2-mu);

                if imag(ttc(1)) == 0
                    ttc = min(ttc);
                    out = stepcounter + ttc;
                    break
                end
            end

            stepcounter = stepcounter + 1;
            A0 = A1;
            B0 = B1;

            if plot_setting == 1
            xlim([-6,6])
            ylim([-4,4])
            plot([r*cos(linspace(0,2*pi,50))+real(A1)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A1)],'green')
            title("simulated walks, history")
            hold on
            plot([r*cos(linspace(0,2*pi,50))+real(B1)]+1i*[r*sin(linspace(0,2*pi,50))+imag(B1)],'blue')
            xlim([-6,6])
            ylim([-4,4])
            hold off
            pause(2^-5)
            end
                    
    end
    
    %if no collision has occured, set:
    if abs( imag(ttc(1)) ) > 0
        out = Inf;
    end
end
end