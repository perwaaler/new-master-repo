function out = find_projection(line1,point)
% finds the point proj_l on line that lies closest to point. line =
% [int,k]. expects complex argument.

% point = -.5 + 2.5i;;
% a = real(point);
% b = imag(point);
% alfa = 1.2;
% beta = 1;
% line1 = [beta,alfa];
% plot(point,'*')
% hold on 
% drawline(line1)

if isreal(point)
    x0=point(1);
    y0=point(2);
else
    x0 = real(point);
    y0 = imag(point);
end

% intercept line1
int1 = line1(1);
% slope line1
k1 = line1(2);

% find parameters of perpendicular line
k2 = -1/k1;
int2 = y0 - k2*x0;

% compute coordinates of projection onto line
x = (int2-int1)/(k1-k2);
y = k2*x + int2;

out = [x,y];

% hold on
% drawline([int2,k2])
% hold on
% plot(x,y,'*')
% xlim([-4,4])
% ylim([-4,4])
% 
% 1/(1+alfa^2)^2*((a*alfa^2-alfa*(b+beta))^2 + (b-beta-alfa*a)^2)
% norm([x,y] - [a,b])^2
% (alfa^-1*a+b-beta)/(alfa^-1+alfa)-x
% (a + alfa*b+  alfa^-1*beta)/(alfa^-1+alfa)-y
% (alfa^-1 + alfa)^-2*( (-alfa*a + b - beta)^2 + (a - alfa^-1*b + alfa^-1*beta)^2)
% 


end

