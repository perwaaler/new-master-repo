
x = linspace(0,0.6*10,50);
y = linspace(0,10,50);
[X,Y] = meshgrid(x,y);

a_t =@(t) logistic_fcn(t, 4,1,[2,5]);
f_d = nan(50,50);
for i=1:50
    for j=1:50
    f_d(i,j) = logistic_fcn(X(i,j), a_t(Y(i,j)), 0.5);
    end
end
f_t = logistic_fcn(Y,     5,  1);


Z = f_d.*f_t;

surf(X,Y,Z)
xlabel('d_{min}')
ylabel('t_{min}')
zlabel('break intensity')

