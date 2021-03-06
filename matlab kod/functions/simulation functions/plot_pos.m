function plot_pos(S, plots, time_diff, decision) 
% plots position of each disc based on their positions in state S.

% unpack variables
A = S(1).pos;
B = S(2).pos;
reactA = S(1).EA;
reactB = S(2).EA;
r = S(1).RUprop.r;
init_x = plots.xinit;
col = plots.col;

% Tadv = -1*-1^time_diff.RU1*time_diff.Tadv;

if nargin == 2
    title_str = sprintf("reactA = %d, reactB = %d", reactA, reactB);
    
elseif nargin == 3
    if isempty(time_diff)
        Tadv = nan;
        TTPC = nan;
    else
        Tadv = time_diff.Tadv;
        TTPC = time_diff.TTPC;
    end
    
    title_str = sprintf("reactA = %d, reactB = %d, Tadv=%d, TTPC=%d", ...
        reactA, reactB, Tadv, TTPC);
    
elseif nargin == 4
    if isempty(time_diff)
        Tadv = nan;
        TTPC = nan;
    else
        Tadv = time_diff.Tadv;
        TTPC = time_diff.TTPC;
    end
    
    col = decision2color(decision,plots);
    title_str = sprintf("reactA = %d, reactB = %d, TTPC=%.0f, Tadv=%.1f, decision=%s", ...
        reactA, reactB, TTPC, Tadv, dec2str(decision(1)));
end

xlimit = [-init_x, init_x];
ylimit = [-7,7];


if plots.enc == 1
    
    xlim(xlimit)
    ylim(ylimit)

    plot(r*exp(1i*linspace(0,2*pi,50)) + A, 'color', col(1))
    
    if plots.EAmode==0
        title(title_str)
    end
    
    hold on
    plot(r*exp(1i*linspace(0,2*pi,50)) + B,'color', col(2))

    xlim(xlimit)
    ylim(ylimit)
    hold off
    
    if max(S.EA)==0
        pause(plots.pause)
    else
        pause(plots.pred_pause)
    end
end
end