function decision = make_decision(time_diff,RUprop,driver_id)

% extract nessecary variables
Tadv = (-1)^(driver_id+1) * time_diff.Tadv;
TTPC = time_diff.TTPC;
aggr = RUprop.aggression(driver_id);

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
[~,decision] = max(z,[],3);

end 

