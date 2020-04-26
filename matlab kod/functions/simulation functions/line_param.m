function line = line_param(A0, Ap)
x0 = real(A0);
y0 = imag(A0);
xp = real(Ap);
yp = imag(Ap);


k = (yp - y0)/(xp - x0);
m = yp - k*xp;
line = [k,m];
end