function color = decision2color(decision,plots)
% function that maps decision(i) (integer) to a color(i)

color(1) = plots.dec_cols(decision(1)+1);
color(2) = plots.dec_cols(decision(2)+1);

end