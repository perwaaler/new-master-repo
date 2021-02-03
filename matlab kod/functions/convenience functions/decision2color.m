function color = decision2color(decision)
% function that maps decision(i) (integer) to a color(i)
color = ["black","black"];
for i=1:2
    
    if decision(i) == 1
        color(i) = "green";
    elseif decision(i) == 2
        color(i) = "cyan";
    elseif decision(i) == 3
        color(i) = "red";
    elseif decision(i) == 4
        color(i) = "yellow";
    elseif decision(i) == 5
        color(i) = "magenta";
    end
    
end

end