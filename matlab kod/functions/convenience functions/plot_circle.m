function a = plot_circle(A,r,col,linewidth)
% function that plots a circle of radius r centered around A

if nargin==2
    hold on
    a = plot(r*exp(1i*linspace(0,2*pi,50)) + A);
    plot(A,'.')
    xlim([-6,6]);
    ylim([-6,6]);
elseif nargin==3
    hold on
    a = plot(r*exp(1i*linspace(0,2*pi,50)) + A,"color",col);
    plot(A,'.')
    xlim([-6,6]);
    ylim([-6,6]);
else
    hold on
    if isempty(col)
        col = 'blue';
    end
    a = plot(r*exp(1i*linspace(0,2*pi,50)) + A,"color",col);
    a.LineWidth = linewidth;
    plot(A,'.')
    xlim([-6,6]);
    ylim([-6,6]);
end

end