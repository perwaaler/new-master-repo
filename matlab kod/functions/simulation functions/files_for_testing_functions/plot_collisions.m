plots.enc = 1;
plots.pred_pause = 0.2;
plots.pause=0.0;
sum(save_enc.type==-2)
ind = find(save_enc.type==-2);

for i=1:length(ind)
    i = ind(i);
%     i = ind(3);
    pause(1)

    for j=1:150
        if isempty(save_enc.states{i}{j})
            break
        end
        Z = save_enc.states{i}{j};
        time_diff = save_enc.time_diff{i}{j};
        decision  = save_enc.decisions{i}(:,j);

        plot_pos(Z,plots,time_diff,decision)
        j=j+1
    end
        
end