
% example of spline interpolation. the y-points contains two more points
% than the x-vector, so spline takes the first and last columns of y as
% endpoint slopes
x = pi*(0:.5:2); 
y = [0  1  0 -1  0  1  0; 
     1  0  1  0 -1  0  1];
pp = spline(x,y);
yy = ppval(pp, linspace(0,2*pi,101));
plot(yy(1,:),yy(2,:),'-b',y(1,2:5),y(2,2:5),'or')
axis equal

% same example, but now we have now condition on the endpoints
y = [ 1  0 -1  0  1 ; 
      0  1  0 -1  0 ];
pp = spline(x,y);
yy = ppval(pp, linspace(0,2*pi,101));
plot(yy(1,:),yy(2,:),'-b',y(1,2:5),y(2,2:5),'or')
axis equal

%%
% spline interpolation using generated path:
state.pos = 1 + 2*1i;
state.theta = 0.1;
state.dtheta = deg2rad(5);
state.speed = 0.3;
state.dspeed = 0.003;
k = 3;

path = gener_pred_path_iter(state,k);
clf
plot_circle(1+2i,0.3,"green")
for i=1:k
plot_circle(path(i),0.3)
end

% u = state.pos - path(1);
% v = path(end) - path(end-1);
% line1 = [real(u);imag(u)];
% line2 = [real(v);imag(v)];
% 
% slope_u = imag(u)/real(u);
% slope_v = imag(v)/real(v);
% 
% x = angle(path); 
% y = zeros(2,k+1);
% theta = angle(path)';
% r = abs(path)';
% y(1,:) = real(path); 
% y(2,:) = imag(path); 
% y = [line1,y,line2];
% y = [[0;.1],y,[0;0]]


theta = angle(path)';
r = abs(path)';
[x,y] = pol2cart(theta,r);

% interpolate with parameter t
t = (1:k);
v = [x,y];
% points to evaluate function

tt = linspace(0,4*k,100);
X = interp1(t,v,tt,'spline');
% plot

pp = spline(t,v');
V = ppval(pp, linspace(1,k*2,100));
plot(V(1,:),V(2,:),'-')
hold on
plot_circle(state.pos,.3)

f_spline =@(x) ppval(pp,x) 

%% Interpolation using interpl

% initial values
state.pos = 1 + 2*1i;
state.theta = 0.1;
state.dtheta = deg2rad(6);
state.speed = 0.1;
state.dspeed = +0.01;

k = 10;
% clf


% generate path
path = gener_pred_path_iter(state,k);
path = path(1:end);

n = length(path);

% plot positions
for i=1:length(path)
plot_circle(path(i),0.3)
end

theta = angle(path)';
r = abs(path)';
[x,y] = pol2cart(theta,r);

% interpolate with parameter t
t = (1:n);
v = [x,y];
% points to evaluate function

tt = linspace(0,4*n,100);
X = interp1(t,v,tt,'spline');
% plot


plot(x,y,'o');
hold on
plot(X(:,1),X(:,2));

for i=1:100
plot_circle(X(i,1) + 1i*X(i,2),0.3,"black")
title(sprintf("time=%0.01f",tt(i)))
pause(0.1)
end

% for t=1:100-n
%     A = X(t+n,1) + 1i*X(t+n,2);
%     plot_circle(A,0.3,"green")
% end



% plot_circle(X(100,1) + 1i*X(100,2),.3,"green")

%% Generating paths and doing trajectory and path-analysis
n = 50;

time.genp  = zeros(1,n);
time.ppfit = zeros(1,n);
time.minim1 = zeros(1,n);
time.minim2 = zeros(1,n);
time.minim2 = zeros(1,n);
time.minim2 = zeros(1,n);
time.pathanal = zeros(1,n);
time.fcn = zeros(1,n);
time.sep = zeros(1,n);

% generate 2 paths:
stateA.pos    = normrnd(-2,1) + normrnd(-2.5,2)*1i;
stateA.theta  = normrnd(deg2rad(10),deg2rad(7));
stateA.dtheta = normrnd(deg2rad(10), deg2rad(6));
stateA.speed  = max([0,gam_rnd(0.3,0.02)-0.05]);
stateA.dspeed = normrnd(0,0.03);  
stateA.RUprop = RUprop;
stateA.id     = 1;

stateB.pos    = normrnd(2,1) + normrnd(1,1)*1i;
stateB.theta  = normrnd(pi,deg2rad(10));
stateB.dtheta = normrnd(0,deg2rad(3));
stateB.speed  = max([0,gam_rnd(0.3,0.02)-0.1]);
stateB.dspeed = gam_rnd(0.01,0.001);
stateB.RUprop = RUprop;
stateB.id     = 2;
plotit        = 1;
pauseL        = 2;

for j=1:n

% stateA = enc.states{stat.frame}(1);
% stateB = enc.states{stat.frame}(2);
states = gen_random_states();
stateA = states{1};
stateB = states{2};

tic
temporal_sep(enc,RUprop,stat,plots);
time.sep(j) = toc;

tic
% path_analysis(ind,stateA,stateB,opt,r,0);
time.fcn(j) = toc;


% k is the total path length. k is number of steps taken, or how far into
% the future the path is generated.
tic

tic
k   = 80;
% ind = floor((1:k) + ((1:k)/4.5).^2 + ((1:k)/9).^3);
% ind = ind(ind<k);
% path2.pos = path2.pos(ind);
paths = gener_pred_path_iter_test(stateA,stateB, k, sparse_ind);
pathA = paths{1};
pathB = paths{2};

% pathA  = gener_pred_path_iter(stateA, k,ind);
% pathB  = gener_pred_path_iter(stateB, k,ind);
k1     = pathA.length;
k2     = pathB.length;
tstop1 = pathA.tstop;
tstop2 = pathB.tstop;
time.genp(j) = toc;

% dists = abs(pathA.padded - pathB.padded);
% loc_extr = find_vector_min(dists);

if plotit==1
clf
for i=1:length(pathA.padded)
    a = plot_circle(pathA.padded(i),0.3,"green");
    if i==length(pathA.s_ind)
        a.Color = "yellow";
    end
    a = plot_circle(pathB.padded(i),0.3,"cyan");
    if i==length(pathB.s_ind)
        a.Color = "yellow";
    end
end
% color starting points red
plot_circle(pathA.pos_dyn(1), 0.3, col2trip("grey"),2);
plot_circle(pathB.pos_dyn(1), 0.3, col2trip("grey"),2);
plot_circle(pathA.padded(pathA.t_proxA), 0.3, col2trip("orange"),2);
plot_circle(pathB.padded(pathB.t_proxB), 0.3, col2trip("orange"),2);
end
tic

%%% translate complex coordinates to cartesian coordinates
x1 = real(pathA.pos_dyn);
y1 = imag(pathA.pos_dyn);

x2 = real(pathB.pos_dyn);
y2 = imag(pathB.pos_dyn);

% parameter values (time)
t1 = pathA.s_ind - 1;
t2 = pathB.s_ind - 1;

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
time.ppfit(j) = toc;
% create piecewise fcn which ensures no crazy behviour for RUs that have
% stopped
ppA =@(t) pw_fcn(ppA,t,tstop1);
ppB =@(t) pw_fcn(ppB,t,tstop2);

% prevent algorithm from selecting times outside of the generated path if
% the acceleration is negative
tstop1 = min(k1,tstop1);
tstop2 = min(k2,tstop2);


%%%% functions that determines search range of fminsearch:
% discourage optimum after RU stops
bound_stop = @(t) (1+exp((t(1)-tstop1)/0.001))*...
                  (1+exp((t(2)-tstop2)/0.001));
% fcn that discourages fmin to search for negative minima
bound_noneg = @(t) (1+exp(-t/0.001));

% define distance functions
d_prox_t =@(t) norm(ppA(t(1)) - ppB(t(2)))*...
                bound_noneg(t(1))*bound_noneg(t(2))*...
                bound_stop(t);
d_min_t  =@(t) norm(ppA(t) - ppB(t))*...
                bound_noneg(t);

%%%% optimization %%%%
% find proximity points, minimum sep. points, and corresponding times
opt = optimset('TolX',0.01,'MaxIter',10,'Display','off');
% t0 = min(t0_cand);
% t0_l = max([0,t0 - 6]);
% t0_r = min([k,t0 + 1]);
tic
% [lbnd,rbnd] = search_range_dmin(loc_extr, ind);
tmin0 = pathA.t_min;
lbnd = tmin0 - 1;
ubnd = tmin0 + 1;
[t_min,d_min] = fminbnd(d_min_t, lbnd, rbnd, opt);
time.minin1(j) = toc;

tic
t_prox = [pathA.t_proxA, pathB.t_proxB];
[t_prox,d_prox] = fminsearch(d_prox_t, t_prox, opt);

[t_prox(1), dproxA] = findTTO(d_prox_t, t_prox, r, 1,[k1,k2],opt);
[t_prox(2), dproxB] = findTTO(d_prox_t, t_prox, r, 2,[k1,k2],opt);

% [t_prox(1,:),d_prox(1)] = fminsearch(d_prox_t, [t_min,t_min], opt);
% t_prox = t_prox(1,:);
% [t_prox(2,:),d_prox(2)] = fminsearch(d_prox_t, flip(t_prox(1,:)), opt);
% [d_prox, I_global] = min(d_prox);
% t_prox = t_prox(I_global,:);
time.minim21(j) = toc; 

tic
% if d_prox < 0.6
%     [t_prox(1), dproxA] = findTTO(d_prox_t, t_prox, r, 1,[k1,k2],opt);
%     [t_prox(2), dproxB] = findTTO(d_prox_t, t_prox, r, 2,[k1,k2],opt);
% end
time.minim22(j) = toc;


% for i=1:2
%     for j=1:2
%     [t_prox,d_prox] = fminsearch(d_prox_t,[] , options);
%     pause(1)
%     end
% end
% 
% 
% a=plot_circle(cart2complex(ppA(t_prox(1))),.3,"black")
% a.LineWidth=2
% a=plot_circle(cart2complex(ppB(t_prox(2))),.3,"black")
% a.LineWidth=2

tic
% find tangent lines at proximity and minimum points
h = 0.0001;
% if k1-1 < t_prox(1)
%     slope_proxA = pathA.slope;
% else
%     slope_proxA = differentiate(ppA, t_prox(1),h);
%     slope_proxA = slope_proxA(2)/slope_proxA(1);
% end

slope_proxA = find_slope(pathA,t_prox(1),ppA,h);
slope_proxB = find_slope(pathB,t_prox(2),ppB,h);
slope_minA  = find_slope(pathA,t_min,ppA,h); 
slope_minB  = find_slope(pathB,t_min,ppB,h); 

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
time.pathanal(j) = toc;

%%%% plotting %%%%
% points to evaluate function (for plotting only)
if plotit==1
    
clf
subplot(121)
plot(dists)

subplot(122)
hold on
for i=1:length(ind)
    a = plot_circle(pathA.padded(i),0.3,"green");
    if i==length(pathA.s_ind)
        a.Color = "yellow";
    end
    a = plot_circle(pathB.padded(i),0.3,"cyan");
    if i==length(ind)
        a.Color = "yellow";
    end
end
% color starting points red
a = plot_circle(pathA.pos_dyn(1), 0.3, col2trip("grey"));
a.LineWidth = 2;
b = plot_circle(pathB.pos_dyn(1), 0.3, col2trip("grey"));
b.LineWidth = 2;

tt1 = linspace(0,k1,200);
tt2 = linspace(0,k2,200);
V1 = ppval(splineA, tt1);
V2 = ppval(splineB, tt2);

% plot spline curve fits
a=plot(V1(1,:),V1(2,:),'-','color','green');
a.LineWidth=1;
a=plot(V2(1,:),V2(2,:),'-','color','blue');
a.LineWidth=1;

% plot proximity points
a=plot_circle(Amin,.3,"black");
a.LineWidth=2;
b=plot_circle(Bmin,.3,"black");
b.LineWidth=2;
drawline([Amin, Amin_pert],'black');
drawline([Bmin, Bmin_pert],'black');
title(sprintf("(%s,%s), t_{min} = %0.01f, t_{prox}=(%0.01f,%0.01f)",...
                            orient(1),orient(2),t_min,ta,tb))

xlim([-10,10])

a = plot_circle(Aprox,.3,col2trip("darkgreen"));
a.LineWidth = 2;
b = plot_circle(Bprox,.3,col2trip("darkblue"));
b.LineWidth = 2;

pause(pauseL)
end

time.total(j) = toc;

end


m.gnp   = median(time.genp);    
m.ppfit = median(time.ppfit);
m.minim1 = median(time.minim1);
m.minim21 = median(time.minim21);
m.minim22 = median(time.minim22);
m.anal  = median(time.pathanal);
m.fcn = median(time.fcn);
m.sep = median(time.sep);
m.tot = sum([m.gnp,m.ppfit,m.minim1,m.minim21,m.minim22,m.anal]);
m

clf
subplot(121)
plot(1:5,[m.gnp, m.ppfit, m.minim1, m.minim21, m.minim22],'o','Color','green');
hold on
plot(6, m.tot,'x','Color',col2trip('darkblue'));
plot(7, m.fcn,'x','Color',col2trip('purple'));
plot(8, m.sep,'x','Color',col2trip('red'));
xlim([0,10])

subplot(122)
hold on
for i=1:length(ind)
    a = plot_circle(pathA.padded_pos(i),0.3,"green");
    if i==length(pathA.s_ind)
        a.Color = "yellow";
    end
    a = plot_circle(pathB.padded_pos(i),0.3,"cyan");
    if i==length(ind)
        a.Color = "yellow";
    end
end
% color starting points red
a = plot_circle(pathA.pos_dyn(1), 0.3, col2trip("grey"));
a.LineWidth = 2;
b = plot_circle(pathB.pos_dyn(1), 0.3, col2trip("grey"));
b.LineWidth = 2;

tt1 = linspace(0,k1,200);
tt2 = linspace(0,k2,200);
V1 = ppval(splineA, tt1);
V2 = ppval(splineB, tt2);

% plot spline curve fits
a=plot(V1(1,:),V1(2,:),'-','color','green');
a.LineWidth=1;
a=plot(V2(1,:),V2(2,:),'-','color','blue');
a.LineWidth=1;

% plot proximity points
a=plot_circle(Amin,.3,"black");
a.LineWidth=2;
b=plot_circle(Bmin,.3,"black");
b.LineWidth=2;
drawline([Amin, Amin_pert],'black');
drawline([Bmin, Bmin_pert],'black');
title(sprintf("(%s,%s), t_{min} = %0.01f, t_{prox}=(%0.01f,%0.01f)",...
                            orient(1),orient(2),t_min,ta,tb))

xlim([-10,10])

a = plot_circle(Aprox,.3,col2trip("darkgreen"));
a.LineWidth = 2;
b = plot_circle(Bprox,.3,col2trip("darkblue"));
b.LineWidth = 2;

pause(pauseL)


