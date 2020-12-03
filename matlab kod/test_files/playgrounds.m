reactivityA =20;
reactivityB =20;
x=-1
Aold = -3 + 0.0001i+x;
A0 = -2+1i+x;
B0 = 6-x + 0.000001i;
Bold = 7-x-1i;
%%

%evasive_step_coulombs_law(A0,B0,stepsize,theta,5, 5)
reactivityA =40;
reactivityB =40;
for i=1:10000

stepsize = [.1 .1];
theta = [degtorad(normrnd(45,2)), degtorad(normrnd(-7,5))];
%Aold = -10 + normrnd(-5,.5)*1i %#ok<*NOPTS>
A0 = normrnd(-5,0.6) + normrnd(-4.5,.8)*1i;
Aold = A0 - stepsize(1)*exp(1i*theta(1));
%Bold = 10 + normrnd(2,0.3)*1i
B0 = normrnd(9,0.3) +  normrnd(2,0.3)*1i;
Bold = B0 + stepsize(2)*exp(1i*theta(2));

% A_real = -xinit;
% B_real = xinit;
% A_im = normrnd(-2.5,sigma_var0);
% B_im = normrnd(0,sigma_var0);
% A0 = A_real + 1i*A_im;
% B0 = B_real + 1i*B_im;
% stepsize = min_stepsize*[1 1] + gamrnd(a_init*[1 1],b_init);   % determines average speed of each vehicle
% stepsize = stepsize/10;
% theta = normrnd([0,0], theta_var_init);
% % %%
% reactivityA =5;
% reactivityB =5;
% x=-0
% Aold = -8 + -8i
% A0   = Aold + 2 + 2i
% Bold =  7  + -1i;
% B0   =  Bold -1.5
% 
% 
% %%
    for k=1:200
        stepsize = [norm(A0-Aold),norm(B0-Bold)];
        
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
        
        delta_speedA = norm(Des_moveA)-norm(predA(1)-A0)
        delta_speedB = norm(Des_moveB)-norm(predB(1)-B0)
        stepsize(1) = stepsize(1) + sign(delta_speedA)*min(0.0006, abs(delta_speedA));
        stepsize(2) = stepsize(2) + sign(delta_speedB)*min(0.0006, abs(delta_speedB));
        stepsize(1) = max(stepsize(1), 0.005)
        stepsize(2) = max(stepsize(2), 0.005)
        %delta_s_max = 0.01;
        delta_theta_max = 2/180*pi;

        % mean and standard deviation of delta theta and delta speed
        mu_thetaA = delta_theta_max*(1 - exp(-1*abs(delta_thetaA)));
        mu_thetaB = delta_theta_max*(1 - exp(-1*abs(delta_thetaB)));
        sigma_thetaA  = degtorad(1) + exp(-abs(delta_thetaA)*5);
        sigma_thetaB  = degtorad(1) + exp(-abs(delta_thetaB)*5);
        logn_parA = lognormal_par(mu_thetaA, sigma_thetaA);
        logn_parB = lognormal_par(mu_thetaB, sigma_thetaB);

        % mu_thetaA/(2*pi) * 360
        % mu_thetaB/(2*pi) * 360
        change_thetaA = sign(delta_thetaA)*lognrnd(logn_parA(1), logn_parA(2));
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

        clf
        plot(A1,'*','color','magenta')
        hold on 
        plot(A0,'*','color','red')
        plot(predA,'*','color','blue')

        plot(B0,'o','color','red')
        plot(B1,'o','color','magenta')
        plot(predB,'o','color','blue')

        xlim([-10,10])
        ylim([-10,10])

        Aold=A0;
        Bold=B0;
        A0=A1;
        B0=B1;
        pause(0.001)
    end


end


% step{1} = [A1,B1];
% step{2} = stepsize;
% step{3} = [thetaA_mod, pi - thetaB_mod];
% out = step;