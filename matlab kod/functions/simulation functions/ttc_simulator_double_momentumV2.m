function out = ttc_simulator_double_momentumV2(A0,B0, stepsize, theta, r, driver_prop, plot_setting)
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
    % generate walk
    while real(A0) < real(B0) + 0.9 && imag(A0)<imag(B0)+0.9
                
                % generate new step
                newstep = take_NEA_step(A0, B0, stepsize, theta, driver_prop);
                A1 = newstep{1}(1);    
                B1 = newstep{2}(2);
                stepsize = newstep{2};
                theta = newstep{3};

                % collision has occured:_
                % if real(A1+r)<real(B1-r)
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
                plot([r*cos(linspace(0,2*pi,50))+real(A1)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A1)])
                title("simulated walks, generated")
                hold on
                plot([r*cos(linspace(0,2*pi,50))+real(B1)]+1i*[r*sin(linspace(0,2*pi,50))+imag(B1)])
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