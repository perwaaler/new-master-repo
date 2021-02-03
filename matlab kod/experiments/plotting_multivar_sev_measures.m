for i=1:50
    i;
    p_i = i/10
x = linspace(0,10,100);
y = linspace(0,3,100);
[X,Y] = meshgrid(x,y);

Pcoll = 1/2*exp(-(X/(2)).^p_i - Y/i);%.*exp(-(Y/4));
Z = (Pcoll).^((Y.^p_i+0.01*X));

clf
surf(X,Y,Z)
title(sprintf("p=%g",p_i))
xlabel('T_{adv}')
ylabel('TTR')
zlabel('P(coll)')
pause(0.0   )
end
%% modelling P_collision
Tadv = linspace(0,3,50);
ttr = linspace(0,6,50);
[X,Y] = meshgrid(Tadv,ttr);

Pcoll = 1/2*exp(-(X/1).^2).*exp(-(Y/2));

surf(Tadv,ttr,Pcoll)

xlabel('T_{adv}')
ylabel('TTR')
zlabel('P(fail resolve)')   
%% T_2

Tadv = linspace(0,3,50);
ttr = linspace(0,6,50);
[X,Y] = meshgrid(Tadv,ttr);

Pcoll = exp(-(X/2+Y/3));

surf(Tadv,ttr,Pcoll)

xlabel('T_{adv}')
ylabel('TTR')
zlabel('P(fail resolve)')  



