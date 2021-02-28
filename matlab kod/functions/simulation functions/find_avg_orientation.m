function side = find_avg_orientation(ZZ,RU1,k,t_prev,n_samplePts,RUprop)
% Computes orientation by taking the most common orientation for the path
% as a whole from the point of wiev of the fastest RU. Returns a character
% vector side that contains the side of RUA and RUB in that order.

RU2 = 3 - RU1;
sample_indeces = linspace(t_prev + 1, k, n_samplePts);
sample_indeces = round(sample_indeces);
left_count =  [0,0];
right_count = [0,0];

for i=1:n_samplePts
    ind = sample_indeces(i);
    RU1info = [RU1,ind,ind];
    RU2info = [RU2,ind,ind];
    sideRU1 = find_orientation(RU1info,ZZ,RUprop);
    sideRU2 = find_orientation(RU2info,ZZ,RUprop);
    
    if     sideRU1=="left"
        left_count(1) = left_count(1) + 1;
    elseif sideRU1=="right"
        right_count(1) = right_count(1) + 1;
    end
    
    if     sideRU2=="left"
        left_count(2) = left_count(2) + 1;
    elseif sideRU2=="right"
        right_count(2) = right_count(2) + 1;
    end
    
%     plot_circle(ZZ{ind}(RU1).pos,r,"magenta")
%     plot_circle(ZZ{ind}(RU2).pos,r,"green")
end

if left_count(1)> right_count(1)
    sideRU1 = "left";
else
    sideRU1 = "right";
end

if left_count(2)> right_count(2)
    sideRU2 = "left";
else
    sideRU2 = "right";
end

side(RU1) = sideRU1;
side(RU2) = sideRU2;

end
