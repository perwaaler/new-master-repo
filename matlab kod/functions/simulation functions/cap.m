function x_capped = cap(x,max_absval)
% returns a value maxcapped with the same sign as x, and with absolute
% value less than or equal to the maximum limit value set by max_absval
x_capped = sign(x)*min(max_absval,abs(x));
end