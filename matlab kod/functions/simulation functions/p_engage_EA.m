function p = p_engage_EA(A0,B0,stepsize0,theta0, detec_par)
% function that computes the probability to EA mode to be engaged
shift_d = detec_par(1);
pow_d   = detec_par(2);
shift_t = detec_par(3);
pow_t   = detec_par(4);
amp     = detec_par(5);

v_delta0 = stepsize0(1)*exp(1i*theta0(1)) + stepsize0(2)* exp(1i*theta0(2));
s_delta0 = A0-B0;
norm_v_delta0 = norm(v_delta0);
t_min = -real(s_delta0*conj(v_delta0))/norm_v_delta0;
d_min = norm(s_delta0 + t_min*v_delta0 );
% p is the product of two s-functions.

%p = 1/(1+exp(pow_d*(d_min - shift_d)))* 1/(1+exp(pow_t*(d_min - shift_t)));
p = amp*exp(-pow_d*(d_min - shift_d))*exp(-pow_t*(d_min - shift_t));
p = exp(-pow_t*(norm(A0-B0) - shift_d))
end

