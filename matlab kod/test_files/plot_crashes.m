% script for studying the encounters that leads to crashes. Used to find
% out whether or not the rules for evasive action makes sense

crashIndex = find(encData.type==-2)';
nCrash = length(crashIndex);

crashStates = encData.states(crashIndex);


for ic=crashIndex
    for j=1:encData.length(ic)
        display(ic)
        S = encData.states{ic}{j};
        time_diff = encData.time_diff{ic}{j};
        decision  = encData.decisions{ic}(:,j);
        plot_pos(S, plots, time_diff, decision)
    end
        
end

