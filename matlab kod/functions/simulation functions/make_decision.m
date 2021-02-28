function decision = make_decision(time_diff, RUprop, RU_id)
% function that simulates a attempt at conflict resolution
% based on Tadv and TTPC.

% extract nessecary variables
RU1  = time_diff.RU1;
Tadv = (-1)^(RU_id+1) * time_diff.Tadv;% negated when driver_id==2
TTPC = time_diff.TTPC;
side = time_diff.side;

if side=="left"
    orientation = 1 ;
else
    orientation = 2;
end
decision_angle = (-1)^orientation;

aggr = RUprop.aggression(RU_id);

%%% define alpha %%%
exp_alpha = binopdf([0,1,2],2,aggr^(exp(-0.4*Tadv)));
% modify the "do_nothing" probability.
exp_alpha(2) = exp_alpha(2)*1/(1 + exp(-(TTPC-15)/4));
% set relative probability to engage in emergency behaviour
exp_alpha(4) = 700/(1+exp((TTPC/3 + abs(Tadv))/1.5))*1/(1+exp((TTPC-10)/0.4 ));
% set alfa
alpha = log(exp_alpha);

% set initial decision
z(1,1,:) = [0,1,0,0];

% simulate sample
z = mrf_sim(z,0,alpha,0,4);
[~,decision_type] = max(z,[],3);

decision = [decision_type, decision_angle];
% if decision==4
%     % RU has determined that the situation is getting severe, and takes
%     % more dramatic action
%     if abs(Tadv)<inf
%         % Tadv is finite, meaning there was prediction-path overlap
%         if Tadv<0
%             % attempt to increase the others time-advantage
%             decision=4;
%         else
%             % attempt to increase own time-advantage
%             decision=5;
%         end
%     else
%         % Tadv is infinite, meaning prediction paths didnt overlap
%         
%             
%     end
% 
% end
            

end 

