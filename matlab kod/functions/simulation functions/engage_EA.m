function out = engage_EA(S,margin)
% criteria for entering Evasive Action

if nargin == 1
    margin = [0.3,0.3];
end

A = S(1).pos;
B = S(2).pos;

% conditions for engaging Evasive Action mode
int_state = max(S(1).EA,S(2).EA);
% relative position conditions
ydiff       = imag(A) < imag(B) + margin(1);
xdiff       = real(A) < real(B) + margin(2);


out = int_state * ydiff * xdiff;
end