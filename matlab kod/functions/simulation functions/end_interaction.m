function interaction = end_interaction(S,max_delta,plots)
% function that decides whether interaction is sufficiently safe to end.
% if RUs can follow their default driving behaviour without getting close,
% interaction is ended

tik = 0; % increases by one for each time distance increases.
D0 = norm(S(1).pos - S(2).pos);
D_min = D0;

RUprop = S.RUprop;


while tik < 3

    E_speedA = normrnd(E_speed_A(S, RUprop), RUprop.Espeed_std(1));
    E_thetaA = normrnd(E_theta_A(S, RUprop), RUprop.Etheta_std(1));
    E_speedB = normrnd(RUprop.avg_speed(2),  RUprop.Espeed_std(2));
    E_thetaB = normrnd(pi,                   RUprop.Etheta_std(2));

    % take step and update states
    S(1) = take_step(S(1), E_speedA, E_thetaA, max_delta(1)); 
    S(2) = take_step(S(2), E_speedB, E_thetaB, max_delta(2));
    
    D = norm(S(1).pos - S(2).pos);
    D_min = min(D0,D);
    
%     plot_pos(S,plots)
    
    % if RUs increase their distance 3 consecutive times, algorithm
    % concludes that 
    if D<D0
        tik = 0;
    else
        tik = tik+1;
    end
    
    D0 = D;
       
end

if D_min > 2
    interaction = 1;
else
    interaction = 0;
end

end