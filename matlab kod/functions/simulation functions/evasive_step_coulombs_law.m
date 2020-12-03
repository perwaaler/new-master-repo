function out = evasive_step_coulombs_law(A0,B0,stepsize,theta, reactivityA, reactivityB)

        Aold = A0 - stepsize(1)*exp(1i*theta(1));
        Bold = B0 + stepsize(2)*exp(1i*theta(2));
        velA = A0-Aold;
        velB = B0-Bold;
        predA = A0 + [1 2]*velA;
        predB = B0 + [1 2]*velB;
        rA = abs(predA(1) - predB);
        rB = abs(predB(1) - predA);

        thetaA = angle(velA);
        thetaB = angle(velB);
        % radtodeg(2)

        % calculate charge of vehicles
        QA = abs(velA);
        QB = abs(velB);

        % calculate the force acting on each car
        FA = reactivityA*QA*QB./(rA-0.2).^2 .* (predA(1)-predB);
        FB = reactivityB*QA*QB./(rB-0.2).^2 .* (predB(1)-predA);

        % calucate desired future positions
        Adesired = predA(1) + sum(FA);
        Bdesired = predB(1) + sum(FB);

        % difference between where they are and where they want to be
        Des_moveA = Adesired - A0;
        Des_moveB = Bdesired - B0;

        % calculate angle of desired move
        delta_thetaA = angle(Des_moveA) - angle(predA(1)-A0);
        delta_thetaB = angle(Des_moveB) - angle(predB(1)-B0);

        %delta_s_max = 0.01;
        delta_theta_max = 2/180*pi;
        
        % change speed
        delta_speedA = norm(Des_moveA)-norm(predA(1)-A0);
        delta_speedB = norm(Des_moveB)-norm(predB(1)-B0);
        stepsize(1) = stepsize(1) + sign(delta_speedA)*min(0.02, abs(delta_speedA));
        stepsize(2) = stepsize(2) + sign(delta_speedB)*min(0.02, abs(delta_speedB));
        stepsize(1) = max(stepsize(1), 0.005);
        stepsize(2) = max(stepsize(2), 0.005);

        % mean and standard deviation of delta theta and delta speed
        mu_thetaA = delta_theta_max*(1 - exp(-1*abs(delta_thetaA)));
        mu_thetaB = delta_theta_max*(1 - exp(-1*abs(delta_thetaB)));
        sigma_thetaA  = degtorad(1) + exp(-abs(delta_thetaA)*30);
        sigma_thetaB  = degtorad(1) + exp(-abs(delta_thetaB)*30);
        logn_parA = lognormal_par(mu_thetaA, sigma_thetaA);
        logn_parB = lognormal_par(mu_thetaB, sigma_thetaB);

        % mu_thetaA/(2*pi) * 360
        % mu_thetaB/(2*pi) * 360
%        change_thetaA = sign(delta_thetaA)*lognrnd(logn_parA(1), logn_parA(2));
%         if abs(change_thetaA)>radtodeg(30)
%             break
%         end

        thetaA_mod = thetaA + sign(delta_thetaA)*lognrnd(logn_parA(1), logn_parA(2));
        thetaB_mod = thetaB + sign(delta_thetaA)*lognrnd(logn_parB(1), logn_parB(2));
        % 0.44/(2*pi) * 360
        
        % update
        theta = [thetaA_mod,thetaB_mod-pi];
        A1 = A0 + stepsize(1)*exp(1i*thetaA_mod);
        B1 = B0 + stepsize(2)*exp(1i*thetaB_mod);


step{1} = [A1,B1];
step{2} = stepsize;
step{3} = theta;
out = step;
end

% angle(velA)/(2*pi)*360
% thetaA/(2*pi)*360
% stepA  = stepA + normrnd(delta_theta_max*(1 - exp(-1*abs(delta_thetaA)),  var_step); 
% stepsize(2) = stepsize(2) + normrnd(0,var_step);




