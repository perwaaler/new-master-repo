
x = linspace(0,2,100);
y = linspace(0,20,70);
[X,Y] = meshgrid(x,y);



Z = logistic_fcn(X, 0.2,-0.1, [1,4*0.6]).*logistic_fcn(Y,5,-1,[0.3,1]);

surf(X,Y,Z)
xlabel('max-speed')
ylabel('time-to-divergence')
zlabel('distance threshold')