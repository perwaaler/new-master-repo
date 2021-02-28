function out = complex2cart(x)
% function that takes complex numbers in x and returns a 
% matrix with each column containing the cartesin coordinates
% of each complex number in x.
out(1,:) = real(x);
out(2,:) = imag(x);
end
