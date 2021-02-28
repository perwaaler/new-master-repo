n_sample = 400;
pl = 1  ;
Tadv_vec = linspace(-12,12,30);
TTPC_vec = linspace(-10,15,25);
TTPC_vec = -6
driver_id = 2;

edges = linspace(0.5, 5.5, 6); % Create 20 bins.
% Plot the histogram.


for i=1:length(TTPC_vec)
    time_diff.TTPC = TTPC_vec(i);
    
    for j=1:length(Tadv_vec)
        
        sample = zeros(1,n_sample);
        time_diff.Tadv = Tadv_vec(j);
        
        for k=1:n_sample
            sample(k) = make_decision(time_diff, RUprop, driver_id);
        end
        
        histogram(sample,'Normalization','probability','BinEdges',edges)
        title(sprintf("TTPC=%0.1f, Tadv=%0.1f",time_diff.TTPC,time_diff.Tadv))
        pause(pl)
    end
    
    pause(0.5)
end