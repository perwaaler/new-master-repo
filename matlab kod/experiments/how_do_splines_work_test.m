
x1 = [0,1,2,2];
y1 = [0,.5,2,3];
slope1 = [1;1];
slope2 = [1;2];
v1 = [x1;y1];
v1 = [slope1,v1,slope2];
v1(:,1) = v1(:,1)*1;
v2(:,1) = v2(:,1)*2;

t = 0:3;
pp = spline(t,v1);

x_eval = linspace(-1,5,1000);
spline_eval = ppval(pp, x_eval);
clf
hold on
p = plot([0,1,2],[0,.5,2],'*');
plot(spline_eval(1,:),spline_eval(2,:))
xlim([-1,3])
ylim([-1,3])
p.LineWidth = 1;

%%

tic
% syms x
% y = piecewise(x<-1, 1/x, x>=-1, sin(x)/x);
% f =@(x) subs(y,x);
fminsearch(@(x) pw_fcn(x),3)
toc

%% 

x1 = [0,1,1.3];
y1 = [0,1,1.7];
v1 = [x1;y1];
t = 0:2;
pp = spline(t,v1);

x_eval = linspace(-1,6,1000);

ppfit =@(x) ppval(pp, x);

spline_eval = ppfit(x_eval);
spline_extrap = ppfit([3,4,5]);
% spline_points
clf
hold on
p = plot(x1,y1,'o');
plot(spline_eval(1,:),spline_eval(2,:))
plot(spline_extrap(1,:),spline_extrap(2,:),'*')
xlim([-1,7])
ylim([-1,7])
p.LineWidth = 1;




