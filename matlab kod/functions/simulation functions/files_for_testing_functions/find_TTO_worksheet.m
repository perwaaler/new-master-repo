[t, d] = findTTO(d_prox_t, t_prox, r, Id, k, contact, opt)
[t_prox(1), dproxA] = findTTO(d_prox_t, paths.t_prox, r, 1,[k1,k2],paths.contact, opt);
% computes the distance and time to reach the closest point in pathB
%%
contact = 0;
opt = optimset('TolX',0.01,'MaxIter',5,'Display','iter');
g =@(X) d_prox_t(X);

% search range for t_proxA:
a  = 0;
b  = paths.t_prox(1);
t_prox = paths.t_prox;

% search limits for ball in path B that A first touches:
bl = max([0, paths.t_prox(2)-7]);
bu = min([k2, paths.t_prox(2)+7]);

if contact==1
    stay_close =@(x) exp( d_prox_t([x, t_prox(2)]*2 ));
else
    stay_close =@(x) 1;
end

D  = 2*r;

[t, d] = fminbnd(@(x) ...
         abs( g([x, fminbnd(@(y) g([x,y]),bl,bu,opt)]) - D ),...
                                                 a,b,opt);
                                             
% g([t, fminsearch(@(y) g([t,y]),t_prox(2),opt)]) - D

plot_circle(cart2complex(ppA(t)),.3,'magenta',2);
% plot_circle(cart2complex(ppB(0.0125)),.3,'magenta',2);



                                             
 %%                                            

    
    g =@(X) d_prox_t(X);

    % search limits for tb:
    a  = 0;
    b  = t_prox(2);
    % search limits for ball in pathA that B first touches:
    al = max([0, t_prox(1)-10]);
    au = min([k(1), t_prox(1)+10]);

    D  = 2*r;

    [t, d] = fminbnd(@(y) ...
                    abs( g([fminbnd(@(x) g([x,y]), al,au,opt),y]) -D),...
                                                  a, b, opt);

