function drawline(line_par,col)
% Draws the line defined by the parameters in line_par.
% The parameters in line_par may be:
% 1) y-intercept and slope                                          [int,k]
% 2) point A line should pass through (complex repr.) and slope:      [A,k]
% 3) two complex numbers [A,B] such that line passes through A and B: [A,B]

if imag(line_par(1))~=0 && imag(line_par(2))~=0
    A0 = line_par(1);
    A1 = line_par(2);
    line_par = find_line_par(A0,A1);
elseif imag(line_par(1))~=1 && imag(line_par(2))==0
    % x0*k + b = y0 --> b = y0 - x0*k
    A = line_par(1);
    k = line_par(2);
    line_par = [imag(A) - real(A)*k, k];
end

int = line_par(1);
k = line_par(2);
xx = linspace(-6,6);

if nargin==1
    hold on
    plot(xx,int + k*xx)
else
    hold on
    plot(xx,int + k*xx,"color",col)
end
end