function [end_criter_met,tick] = stop_gen_path(dist,contact,...
             first_zero_step,done_intersecting,tick,end_criter_met,analyze)

%%% EXIT 1: ARE THEY DIVERGING WITHOUT MAKING CONTACT?
if dist.prox(1)<=dist.prox(2)

    if not(contact)
        % paths have not gotten closer
        tick.div = tick.div + 1;
        if tick.div == 5
            end_criter_met = true;
            if analyze
                disp('divergence with no contact')
            end
        end
    end
    
else
        % reset divergence ticker everytime t_prox improves
        tick.div = 0;
end

%%% EXIT 2: HAVE BOTH STOPPED?
if min(first_zero_step)==1
    end_criter_met = true;
    if analyze
        disp('both stopped')
    end
end

%%% EXIT 3: IF THEY'VE INTERSECTED, ARE THEY DONE INTERSECTING?
if done_intersecting.A && done_intersecting.B
    end_criter_met = true;
    if analyze
        disp("finnished intersecting")
    end
end

end