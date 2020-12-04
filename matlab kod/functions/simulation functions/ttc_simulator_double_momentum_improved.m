function out = ttc_simulator_double_momentum_improved(A0,B0,stepsize0,theta0,min_stepsize,sigma_step, r, thetavar,stability_fac, plot_setting)
% out = ttc_simulator(A0,B0,stepsize0,min_stepsize,sigma_step, r, thetavar)
% Generates a random walk starting at position A0 and
% B0 respectively, and then computes Time To Collision. TTC is set to inf
% if collision never occurs.
% A0 and B0 given current positions, and theta describe current direction.

stepcounter = 0;
ttc = 1i;
D = norm(A0-B0);

% if they have already collided, run:
if D < 2*r
    % get pre-collision positions based on arguments
    B1 = B0 + stepsize0(2)*exp(1i*theta0(2));
    A1 = A0 - stepsize0(1)*exp(1i*theta0(1));
    Adiff = A1-A0;
    Bdiff = B1-B0;
    eta = real(conj(Adiff - Bdiff)*(A0-B0))/norm(Adiff-Bdiff)^2;
    mu = (norm(A0 - B0)^2 - 4*r^2)/norm(Adiff - Bdiff)^2;
    ttc = -eta + [-1 1]*sqrt(eta^2-mu);
    ttc = -max(ttc); % ttc will now be negative, since time must be reversed to find time of collision
    out = ttc;
else
    % generate walk
    while real(A0) < real(B0)
                
                % generate new step
                theta0 = theta0 + normrnd(-stability_fac*theta0,thetavar);
                stepsize0 = stepsize0 + normrnd([0,0],sigma_step);
                stepsize0(1) = max(min_stepsize, stepsize0(1));
                stepsize0(2) = max(min_stepsize, stepsize0(2));
                A1 = A0 + stepsize0(1)*exp(1i*theta0(1));
                B1 = B0 - stepsize0(2)*exp(1i*theta0(2));        

                % collision has occured
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
                
%                 if real(A1) < real(B1)
%                     if 
%                 end
                
                if plot_setting == 1
                xlim([-6,6])
                ylim([-4,4])
                plot([r*cos(linspace(0,2*pi,50))+real(A1)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A1)])
                title("simulated walks")
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