function out = sevme_str(sev_ind)
% function that translates the transformation index to a string that can be
% inserted into plot title and legend.
sev_measure = cell(1,3);
sev_measure{1} = 'stoch. TTC';
sev_measure{2} = 'TTC';
sev_measure{3} = 'min. dist.';

out = sev_measure{sev_ind};
end