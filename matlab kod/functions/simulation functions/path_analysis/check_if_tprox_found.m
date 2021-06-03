function [tprox_found,tick,dist] = check_if_tprox_found(tick,dist)                                  
% returns 1 if stopping criteria for path generation are reached.
% path generation stops if the paths are categorized as diverging without
% ever having touched OR if they are done overlapping.

tprox_found = false;

%%% has approximation of t_prox been found? %%%
if dist.prox(2)==dist.prox(1)
    tick.t_prox = tick.t_prox + 1;
    % the drivers have made contact in the past, and are no longer in
    % contact
    if tick.t_prox >= 5
        tprox_found = true;
    end
end

end