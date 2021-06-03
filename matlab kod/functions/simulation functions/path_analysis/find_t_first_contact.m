function t_low_bnd = find_t_first_contact(d_prox_t,t_prox,r)
ta = t_prox(1);
tb = t_prox(2);
d_thr = 2*r;
opt = optimset('maxIter',3,'Display','off');


d = 0;
j = 0;
% distance from A(ta) to point in B closest to A(ta)
while d<d_thr && ta-j>0
    j = j + 1;
    d = d_prox_t([ta-j, fminsearch(@(y) d_prox_t([ta-j,y]),tb,opt)]);

    if ta-j<0
        break
    end
end
t_low_bnd(1) = max(0,ta - j);
    
d = 0;
j = 0;
% distance from A(ta) to point in B closest to A(ta)
while d<d_thr && tb-j>0
    j = j + 1;
    d = d_prox_t([ fminsearch(@(x) d_prox_t([x,tb-j]),ta,opt),tb-j]);
    if tb-j<0
        break
    end
end
t_low_bnd(2) = max(0,tb - j);
    

end