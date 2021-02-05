function color = decision2color(decision,plots)
% function that maps decision(i) (integer) to a color(i)

color(1) = plots.dec_cols(decision(1)+1);
color(2) = plots.dec_cols(decision(2)+1);


% for i=1:2
%     if        decision ==0
%         color(i) = "black";
%     elseif decision(i) == 1
%         color(i) = "green";
%     elseif decision(i) == 2
%         color(i) = "cyan";
%     elseif decision(i) == 3
%         color(i) = "red";
%     elseif decision(i) == 4
%         color(i) = "yellow";
%     elseif decision(i) == 5
%         color(i) = "magenta";
%     end
%     
% end

end