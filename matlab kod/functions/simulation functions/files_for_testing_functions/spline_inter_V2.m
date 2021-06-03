

plotit        = 1;
pauseL        = 0;

n=1000;
t.sep = zeros(1,n);
t.tot = zeros(1,n);
t.gen  = zeros(1,n);
t.slope = zeros(1,n);
t.opt1 = zeros(1,n);
t.opt2 = zeros(1,n);
% stateA = enc.states{stat.frame}(1);
% stateB = enc.states{stat.frame}(2);
% states = gen_random_states();
% stateA = states{1};
% stateB = states{2};
% for j=1:30
sparse_ind = floor((1:k) + ((1:k)/3.0).^2 + ((1:k)/7).^3);
% ind = ind(ind<k);
for j=1:n
% k is the total path length. k is number of steps taken, or how far into
% the future the path is generated.
k   = 80;

states = gen_random_states();
stateA = states{1};
stateB = states{2};


enc.states{1} = stateA;
enc.states{2} = stateB;
enc.frame = 1;
RUprop = stateA.RUprop;


tic
paths = gener_pred_path_iter_test(stateA,stateB,k,sparse_ind,false,false);
pathA = paths.A;
pathB = paths.B;
t.gen(j) = toc;

k1     = paths.n_rel(1);
k2     = paths.n_rel(2);

if plotit==1
clf
for i=1:length(pathA.full)
    a = plot_circle(pathA.full(i),0.3,"green");
    if i==length(pathA.s_ind)
        a.Color = "yellow";
    end
end
for i=1:length(pathB.full)
    a = plot_circle(pathB.full(i),0.3,"cyan");
    if i==length(pathB.s_ind)
        a.Color = "yellow";
    end
end
% color starting points red
plot_circle(pathA.pos_sparse(1), 0.3, col2trip("grey"),4);
plot_circle(pathB.pos_sparse(1), 0.3, col2trip("grey"),4);
plot_circle(pathA.full(paths.t_prox(1)+1), 0.3, col2trip("orange"),2);
plot_circle(pathB.full(paths.t_prox(2)+1), 0.3, col2trip("orange"),2);
plot_circle(pathA.full(paths.t_min+1), 0.3, col2trip("orange"),2);
plot_circle(pathB.full(paths.t_min+1), 0.3, col2trip("orange"),2);
end

tic
%%% translate complex coordinates to cartesian coordinates
x1 = real(pathA.pos_sparse);
y1 = imag(pathA.pos_sparse);

x2 = real(pathB.pos_sparse);
y2 = imag(pathB.pos_sparse);

% parameter values (time)
t1 = ind2time(pathA.s_ind);
t2 = ind2time(pathB.s_ind);

% corresponding output values (2d-positions)
v1 = [x1;y1];
v2 = [x2;y2];

% special case when the RU is standing still
if size(v1,2)==1
    v1 = [v1,v1]; %#ok<*AGROW>
    t1 = [0,1];
end
if size(v2,2)==1
    v2 = [v2,v2];
    t2 = [0,1];
end


% fit splines
splineA = makima(t1,v1);
splineB = makima(t2,v2);
ppA =@(t) ppval(splineA,t);
ppB =@(t) ppval(splineB,t);
% create piecewise fcn which ensures no crazy behviour for RUs that have
% stopped

t_relA = paths.rel_time(1);
t_relB = paths.rel_time(2);
ppA =@(t) pw_fcn(ppA,t,t_relA);
ppB =@(t) pw_fcn(ppB,t,t_relB);
t.spline(j) = toc;

% prevent algorithm from selecting times outside of the generated path if
% the acceleration is negative

%%%% functions that determines search range of fminsearch:
% discourage optimum after RU stops
bound_stop = @(t) (1+0.2*exp((t(1)-t_relA)/0.01))*...
                  (1+0.2*exp((t(2)-t_relB)/0.01));
% fcn that discourages fmin to search for negative minima
bound_noneg = @(t) (1+0.2*exp(-t/0.01));

% bound_stop = @(t) 1;
% bound_noneg = @(t) 1;


% define distance functions
d_prox_t =@(t) norm(ppA(t(1)) - ppB(t(2)))*...
                bound_noneg(t(1))*bound_noneg(t(2))*...
                bound_stop(t);
d_min_t  =@(t) norm(ppA(t) - ppB(t))*...
                bound_noneg(t);

            

%%%% optimization %%%%
% find proximity points, minimum sep. points, and corresponding times
opt = optimset('TolX',0.01,'MaxIter',5,'Display','off');
% t0 = min(t0_cand);
% t0_l = max([0,t0 - 6]);
% t0_r = min([k,t0 + 1]);
% [lbnd,rbnd] = search_range_dmin(loc_extr, ind);
tic
tmin0 = paths.t_min;
lbnd = max(0,tmin0 - 1);
ubnd = min(  tmin0 + 1,max(t_relA,t_relB));
[t_min,d_min] = fminbnd(d_min_t, lbnd, ubnd, opt);
t.opt1(j) = toc;

% for j=1:30
tic
t_prox=[0,0];
if paths.contact
    [t_prox(1), dproxA] = findTTO(d_prox_t, paths, r, 1,opt);
    [t_prox(2), dproxB] = findTTO(d_prox_t, paths, r, 2,opt);
else
    [t_prox, dproxA] = findTTO(d_prox_t, paths, r, 1,opt);
end
t.opt2(j) = toc;
% end
% median(t.opt2);
    


tic
% find tangent lines at proximity and minimum points
h = 0.0001;

slope_proxA = find_slope(paths,1,t_prox(1),ppA,h);
slope_proxB = find_slope(paths,2,t_prox(2),ppB,h);
slope_minA  = find_slope(paths,1,t_min,ppA,h); 
slope_minB  = find_slope(paths,2,t_min,ppB,h); 

% compute proximity points and minimality points
Aprox = cart2complex(ppA(t_prox(1)));
Bprox = cart2complex(ppB(t_prox(2)));
Amin  = cart2complex(ppA(t_min));
Bmin  = cart2complex(ppB(t_min));

Aprox_pert = Aprox + 1 + 1i*slope_proxA;
Bprox_pert = Bprox + 1 + 1i*slope_proxB;
Amin_pert  = Amin + 1 + 1i*slope_minA;
Bmin_pert  = Bmin + 1 + 1i*slope_minB;

orient(1) = find_side(1 + 1i*slope_minA,Bmin-Amin);
orient(2) = find_side(1 + 1i*slope_minB,Amin-Bmin);
t.slope(j)=toc;


tic
time_sep = temporal_sep(enc,RUprop,stat,plots,...
                                 8, 80, 0);
t.sep(j) = toc;
    
%%%% plotting %%%%
% points to evaluate function (for plotting only)
if plotit==1
    
clf
hold on
for i=1:length(pathA.full)
    a = plot_circle(pathA.full(i),0.3,"green");
    if i==length(pathA.s_ind)
        a.Color = "yellow";
    end
end
for i=1:length(pathB.full)
    a = plot_circle(pathB.full(i),0.3,"cyan");
    if i==length(ind)
        a.Color = "yellow";
    end
end
% color starting points red
plot_circle(pathA.full(1), 0.3, col2trip("grey"),4);
plot_circle(pathB.full(1), 0.3, col2trip("grey"),4);

tt1 = linspace(0,k1,200);
tt2 = linspace(0,k2,200);
V1 = ppval(splineA, tt1);
V2 = ppval(splineB, tt2);

% plot spline curve fits
plot(V1(1,:),V1(2,:),'-','color','green','LineWidth',1);
plot(V2(1,:),V2(2,:),'-','color','blue' ,'LineWidth',1);

% plot minimum points
plot_circle(Amin,.3,"black",3);
plot_circle(Bmin,.3,"black",3);
drawline([Amin, Amin_pert],'black');
drawline([Bmin, Bmin_pert],'black');
title(sprintf("(%s,%s), t_{min} = %0.01f, t_{prox}=(%0.01f,%0.01f)",...
                            orient(1),orient(2),t_min,ta,tb))

xlim([-10,10])

% draw proximity points
plot_circle(Aprox,.3,col2trip("pink"),2);
plot_circle(Bprox,.3,col2trip("red"),2);

end
end


t.tot = t.gen + t.opt2 + t.opt1 + t.spline + t.slope;

t.sep = t.sep(10:end);
t.gen = t.gen(10:end);
t.spline = t.slope(10:end);
t.opt1 = t.opt1(10:end);
t.opt2 = t.opt2(10:end);
t.slope = t.slope(10:end);
t.tot = t.gen + t.opt2 + t.opt1 + t.spline + t.slope;

m.sep=median(t.sep)
m.gen = median(t.gen)
m.opt2 = median(t.opt2)
m.opt1 = median(t.opt1)
m.spline = median(t.spline)
m.slope = median(t.slope)
m.tot = median(t.tot);
avg.sep = mean(t.sep)
avg.gen = mean(t.gen)
avg.opt2 = mean(t.opt2)
avg.opt1 = mean(t.opt1)
avg.spline = mean(t.spline)
avg.slope = mean(t.slope)
avg.tot = mean(t.tot)

clf
plot([1,2,3,4,5,6,7],[m.sep,m.tot,m.gen,m.spline,m.opt1,m.opt2,m.slope],...
    'Color','green')
hold on
plot([1,2,3,4,5,6,7],...
    [avg.sep,avg.tot,avg.gen,avg.spline,avg.opt1,avg.opt2,avg.slope],...
    'Color','blue')

legend('median','mean')

avg.tot/avg.sep

