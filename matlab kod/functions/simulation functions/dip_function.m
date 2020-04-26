function speed = dip_function(A0, d, c, p, l, amp)
% defines a function symmetric around d, which has minimum at d, and goes
% takes maximum as x goes to +/- infty.
speed = (s_function(real(A0), d-c, p) + s_function(real(A0), d+c, -p) + l)*amp;
end