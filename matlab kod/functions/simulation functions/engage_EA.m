function out = engage_EA(S,stat,margin)
% criteria for entering Evasive Action

if nargin == 2
    margin = [0,0.5];
end

A = S(1).pos;
B = S(2).pos;
out = stat.int==1 && real(A) < real(B) && imag(A) < imag(B)+margin(2);

end