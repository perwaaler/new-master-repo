function [succ_rate, fail_rate] = est_succ_fail_rate(matrix1,matrix2)
% function that estimates the rate at which method fails to correctly
% predict which data belongs to the safe trafic location. matrix1 and 
% matrix2 contains probability estimates or return level estimates from 
% situation of safety level 1 and 2 respectively. Returns estimates of 
% failure rates for each threshold.

succ_rate = sum(matrix2 > matrix1)/500;   
fail_rate = sum(matrix2 < matrix1)/500;

end