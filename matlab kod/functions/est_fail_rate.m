function [succ_rate, fail_rate] = est_succ_fail_rate(p1,p2)
% function that estimates the rate at which method fails to correctly
% predict which data belongs to the safe trafic location. p1 and p2
% contains probability estimates from situation of safety level 1 and 2
% respectively. Returns estimates of failure rates for each threshold.

succ_rate = sum(p1 < p2)/500;   
fail_rate = sum(p1 > p2)/500;

end