function out = continue_simulation(S,xinit,min_speed,margin)
% criteria for continuing simulation

A = S(1).pos;
B = S(2).pos;

% set default arguments
if nargin == 2
    margin = [1.0 1.0];
    min_speed = 1e-4;
    
elseif nargin == 3
    margin = [1.0 1.0];
end

% boundary conditions
x_boundaryA = real(A) < xinit;
y_boundaryA = imag(A) < 4;
x_boundaryB = real(B) > -xinit;
% relative position conditions
ydiff       = imag(A) < imag(B) + margin(1);
xdiff       = real(A) < real(B) + margin(2);
% speed condition
speed    = max(S(1).speed,S(2).speed) > min_speed;


% all above must be satisfied to continue
out = x_boundaryA * y_boundaryA * x_boundaryB * ydiff * xdiff * speed;


end