Tadv_thr = 5;
% how large does change in Tadv have to be to indicate communication?
delta_Tadv_thr = 3;
% determines the threshold below which the situation is uncomfertable
TTPC_thr = 10;
% determines propensity for risktaking behaviour
aggression = beta_rnd(0.5,0.15);

%%
decision = [];
clear_Tadv    = Tadv > Tadv_thr;
clear_Tdisadv = Tadv < Tadv_thr;
first_decision = stat.decision(1)==0;
large_TTPC     = TTPC > TTPC_thr;

if first_decision==0
    delta_Tadv = stat.Tadv(2) - stat.Tadv(1);
end

decision_vec = ["let_other_first","do_nothing","go_first",...
                "left_and_break","right_and_break"];


if clear_Tadv==1 
    % clear who has time advantage
    if first_decision==1
        decision='go_first';
        
    elseif delta_Tadv > delta_Tadv_thr
        decision='go_first';
    elseif abs(delta_Tadv) < delta_Tadv_thr
        decision='go_first';
    else
        decision='let_other_pass';
    end

else
    % not clear who has time advantage
    if large_TTPC==1
        % there is still safe amount of time to solve conflict
       
        if first_decision==1
            alpha = log(binopdf([0,1,2],2,aggression));
            z0(1,1,:) = [1,0,0];
            decision = mrf_sim(z0,0,alpha,0,3);
            [~,decision] = max(decision,[],3);
            decision = decision_vec(decision);
            
        elseif first_decision==0
            % small Tadv>0, large decision time, have made decision before
            if delta_Tadv > delta_Tadv_thr
                % other RU is signalling he wants to let you go first
                decision = 'go_first';
            elseif abs(delta_Tadv) < delta_Tadv_thr
                % no clear signal from other driver, make decision
                alpha = log(binopdf([0,1,2],2,aggression));
                z0(1,1,:) = [1,0,0];
                decision = mrf_sim(z0,0,alpha,0,3);
                [~,decision] = max(decision,[],3);
                decision = decision_vec(decision);
                % choose either to go first, let other pass, or do nothing
            elseif delta_Tadv < delta_Tadv_thr 
                % other driver is signalling he wants to go first, so let
                % him pass
                decision='let_other_pass';
            end
            
        end
                  
    else
        % Small Tadv and small TTCP, situation is now severe.    
        decision = 'left_and_break';
    end
            
    
end  
        
decision
