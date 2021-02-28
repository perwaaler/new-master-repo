function line_par = find_line_par(A0,A1)
% function that takes two points in the complex plane and computes the
% y-intercept and slope of the line that pass through them.

x0 = real(A0);
x1 = real(A1);
y0 = imag(A0);
y1 = imag(A1);

slope = (y1-y0)/(x1-x0);
int = y0 - slope*x0;

line_par = [int,slope];
end