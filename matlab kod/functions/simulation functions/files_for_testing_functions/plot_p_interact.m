x = linspace(0,20,50);
y = linspace(0,3,50);
[X,Y] = meshgrid(x,y);


for i=1:20
b_speed = -2*i/10;
beta = [beta_EA.int;beta_EA.dist;beta_EA.speed];
pea = RUprop.alert(1)./(1 + exp(beta_EA.int + beta_EA.dist*X + b_speed*Y));


surf(X,Y,pea)
title(sprintf('b_{speed} = %g',b_speed))
xlabel('distance')
ylabel('speed')
zlabel('p_{react}')
pause(.5)
end
