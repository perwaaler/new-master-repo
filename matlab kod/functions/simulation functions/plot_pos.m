function plot_pos(S, plots, init_x, col, Tadv) 
% plots position of each disc based on their positions in state S.

% unpack variables
A = S(1).pos;
B = S(2).pos;
reactionA = S(1).EA;
reactionB = S(2).EA;
r = S(1).RUprop.r;
% Tadv = -1*-1^time_diff.RU1*time_diff.Tadv;

if nargin == 2
    init_x = 4;
elseif nargin == 3
    col = ["red", "blue"];
    Tadv = nan; 
elseif nargin == 4
    Tadv = nan;
end

xlimit = [-init_x, init_x];
ylimit = [-5,5];

if plots.enc == 1
    
    xlim(xlimit)
    ylim(ylimit)

    plot(r*exp(1i*linspace(0,2*pi,50)) + A, 'color', col(1))
    
    Tadv = num2str(Tadv);
    title(sprintf("reactionA = %d, reactionB = %d, Tadv=%s", reactionA,reactionB,Tadv))
    hold on 

    plot(r*exp(1i*linspace(0,2*pi,50)) + B,'color', col(2))

    xlim(xlimit)
    ylim(ylimit)
    hold off
    pause(plots.pause)
    
    if norm(A - B) < 2*r
        hold on; title("collision"); hold off
    end
end
end