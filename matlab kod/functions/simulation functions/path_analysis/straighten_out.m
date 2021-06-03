function free_theta = straighten_out(stateA,stateB,...
                                    thetaA_0,thetaB_0,free_theta)
% checks if the road users should straighten out in order to prevent them
% from going in circle

%%% make sure neither RU goes in a circle %%%
if abs(stateA.theta - thetaA_0)>=2.6180
    % force theta to be fixed after turning more than 2.6
    free_theta(1) = 0;
elseif abs(stateB.theta - thetaB_0)>=2.6180
    free_theta(2) = 0;
end

end