
% n = 10;
% x0 = linspace(l_bound,u_bound,n);
% y0 = linspace(l_bound,u_bound,n);
% linspace(l_bound,u_bound,2)
% trial = zeros(n,n);
% 
% for i=1:n
%     for j=1:n
%         trial(i,j) = f([x0(i),y0(j)]);
%     end
% end
% [~,I] = min(trial(:),[],'all') 
% 
%     
% [T(i,:),d(1)] = fminsearch(d_prox_t, [3,10], options);

problem = createOptimProblem('fmincon',...
    'objective',@(x)f(x),...
    'x0',[1,1],'options',...
    optimoptions(@fmincon,'Algorithm','sqp','Display','off'));

%%
p(1,2) = .31; p(2,1)=p(1,2);
p(3,4) = .31; p(4,3)=p(3,4);
p(1,3) = .2;  p(3,1)=p(1,3);
p(1,4) = .14; p(4,1)=p(1,4);
p(2,3) = .03; p(3,2)=p(2,3);
p(2,4) = .01; p(4,2)=p(2,4);

pi = zeros(1,4);
for i=1:4
    ind  = 1:4;
    ind(i) = [];
    pi(i) = sum(p(i,ind));
    pause(1)
end

%% testing the constrained minimization function
 %#ok<*ASGLU>

opt1 = optimset('TolX',10^-1,'MaxIter',5);
opt2 = optimset('TolX',10^-1,'MaxIter',30);


n=500;

tnew = 1:n;
for i=1:n
tic
[tmin_new,dmin_new] = fminbnd(d_min_t,0,k,opt1);
tnew(i) = toc;
end
plot(tnew(10:end))
mnew = median(tnew)


told=1:n;
for i=1:n
tic
[tmin_old,dmin_old]   = fminsearch(d_min_t, t_min, opt2);
told(i) = toc;
end
plot(told(10:end))
mold = median(told)


[(tmin_new/13.4179),(tmin_old/13.4179)]
mnew/mold %#ok<*NOPTS>
% the one dimensional optimizator is way better! 70% faster!

%% testing the bounded optimization with several variables

A = [[-1,0];[0,-1];[1,0];[0,1]]
b = [0;0;k1;k2]
[x1,x2] = fmincon(d_prox_t,[1,1],A,b)



options = optimset('TolX',0.1,'MaxIter',200);

n=100;

told=1:n;
for i=1:n
tic
[t_m,d_m] = fminsearch(d_prox_t, [1,1], options);
told(i) = toc;
end
plot(told(10:end))
mold = median(told)


tnew = 1:n;
for i=1:n
tic
% A = [[-1,0];[0,-1];[1,0];[0,1]];
% b = [0;0;30;30];
x0 = [15,15];
A = [];
b = [];
UB = [30;30];
LB = [10;10];
[x1,x2] = fmincon(d_prox_t,[1,1],[],[],[],[],LB,UB,[],options);
tnew(i) = toc;
end
plot(tnew(10:end))
mnew = median(tnew)

mnew/mold %#ok<*NOPTS>
% Conclusion: bounded optimization is better if the search range
% is narrow

%% speed comparison; who can compute ttpath fastest?
t2 = 3
f = @(t1) d_prox_t([t1,t2])
f1 =@(t1) fminbnd(@(t1)    d_prox_t([t1,t2]),0,k, options)
f2 =@(t1) fminsearch(@(t1) d_prox_t([t1,t2]),1,   options)

t = zeros(2,100);
%%%%% 1 %%%%%
for i=1:100
tic
f1(1);
t(1,i) = toc;
end
%%%%% 2 %%%%%
for i=1:100
tic
f2(1);
t(2,i) = toc;
end
median(t(2,:))/median(t(1,:))
% still fminbnd is much faster
%%
options = optimset("TolX", 10^-3)
t_prox
f =@(t1) fminbnd(@(t1)  d_prox_t([t1,15]),0,k, options);
f(20)
%%
% find closest point on other path

f = @(t1) d_prox_t( [fminbnd(@(t1) d_prox_t([t1,t2],0,k)), t2] );

g =@(t2) fminbnd(@(t2) f,...
                    0,k, options);

tic               
fminsearch(g,11,options)
toc
%% finding ta
optx = optimset("TolX", 10^-1, "MaxIter",20)
opty = optimset("TolX", 10^-1, "MaxIter",20)

g =@(X) d_prox_t(X);
a=0;
b=t_prox(1);
bl = max([0, t_prox(2)-20]);
bu = min([k2, t_prox(2)+20]);

tic
[xmin, fcn_min] = fminbnd(@(x) abs(g([x, fminbnd(@(y) g([x,y]), bl,bu,opt)]) -0.6 ),...
                                            a,b,optx);

toc
[t_prox(1), dproxA] = findTTO(d_prox_t, t_prox, r, 1,[k1,k2],opt)

a=plot_circle(cart2complex(ppA(xmin)),r);
a.LineWidth=2;
a.Color=col2trip("darkred");

%% finding tb

opt = optimset("TolX", 10^-1, "MaxIter",100)



g =@(X) d_prox_t(X);
% a=t_prox(2);
% b=k;
tic
a=0;
b  = t_prox(2);
al = max([0, t_prox(1)-20]);
au = min([k, t_prox(1)+20]);

[ymin, fcn_min] = fminbnd(@(y) ...
                    abs(g([fminbnd(@(x) g([x,y]), bl,bu,opt),y]) -0.6 ),...
                                            a,b,opt);

toc
ymin
a=plot_circle(cart2complex(ppB(ymin)),r);
a.LineWidth=2;
a.Color=col2trip("darkred");
%%
stateA = enc.states{stat.frame}(1);
stateB = enc.states{stat.frame}(2);
plotOpts.pause = 1;
plotOpts.PL = 1;
n=50;
opt = optimset('TolX',0.1,'MaxIter',10,'Display','off');

t2sample = zeros(1,n);
for i=1:n
tic
path_analysis(ind,stateA,stateB,opt,r,0);
t2sample(i) = toc;
end
m2 = median(t2sample)


t1sample = zeros(1,n);
for i=1:n
tic
temporal_sep(enc,RUprop,stat,plots);
t1sample(i) = toc;
end
m1 = median(t1sample)

t3sample = zeros(1,n);
for i=1:n
tic
gener_pred_path_iter_test(stateA,stateB, 80, ind);
t3sample(i) = toc;
end
m3 = median(t3sample)

t4sample = zeros(1,n);
for i=1:n
tic
gener_pred_path_iter(stateA, k, ind);
t4sample(i) = toc;
end
m4 = median(t4sample)



m2/m1

%%
n=100;
t.sep = nan(1,n);
t.spline = nan(1,n);
genInfo.plotP = 4;
for j=1:n
k   = 80;
sparse_ind = floor((1:k) + ((1:k)/4).^2 + ((1:k)/7).^3);
% ind = ind(ind<k);
states = gen_random_states();
stateA = states{1};
stateB = states{2};

RUprop = stateA.RUprop;
genInfo.nIter = 12;
tic
spline_pred_path(stateA,stateB,k,sparse_ind,genInfo)
t.spline(j) = toc;
end

for j=1:n
k   = 80;
% ind = floor((1:k) + ((1:k)/4.5).^2 + ((1:k)/9).^3);
% ind = ind(ind<k);
states = gen_random_states();
stateA = states{1};
stateB = states{2};

enc.states{1} = stateA;
enc.states{2} = stateB;
enc.frame = 1;
RUprop = stateA.RUprop;

tic
time_sep = temporal_sep(enc,RUprop,stat,plots,...
                                 8, 80, 0);
t.sep(j) = toc;
end



mean(t.spline(10:end))/mean(t.sep(10:end))

clf
plot(t.spline(10:end))
hold on
plot(t.sep(10:end))


%%
n=100;
t = nan(1,n);


for j=1:n
tic
t_prox = paths.t_prox;
x = normrnd(2,1000);
abs(x);
t(j) = toc;
end
mean(t)

for j=1:n
tic
t_prox=[0,0];
if paths.contact
    [t_prox(1), dproxA] = findTTO(d_prox_t, paths, r, 1,opt);
    [t_prox(2), dproxB] = findTTO(d_prox_t, paths, r, 2,opt);
else
    [t_prox, dproxA] = findTTO(d_prox_t, paths, r, 1,opt);
end
t_opt2(j) = toc;
end

plot(t_opt2(10:end))


for j=1:n
tic
t_prox=[0,0];
if paths.contact
    [t_prox(1), dproxA] = findTTO_V2(d_prox_t, paths, r, 1,opt);
    [t_prox(2), dproxB] = findTTO_V2(d_prox_t, paths, r, 2,opt);
else
    [t_prox, dproxA] = findTTO(d_prox_t, paths, r, 1,opt);
end
t_opt3(j) = toc;
end

mean(t_opt2)/mean(t_opt3)





