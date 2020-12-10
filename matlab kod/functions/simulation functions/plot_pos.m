function plot_pos(S,D, plotter, xinit, r) 
% plots position of each disc based on their positions in state S.

% unpack variables
A = S(1).pos;
B = S(2).pos;
reactionA = S(1).EA;
reactionB = S(2).EA;
clf
if plotter.enc == 1
    
    hold on 
    xlim([-xinit,xinit])
    ylim([-4,4])
    
    plot(r*exp(1i*linspace(0,2*pi,50)) + A)
    title(sprintf("reactionA = %d, reactionB = %d", reactionA,reactionB))
    
    plot(r*exp(1i*linspace(0,2*pi,50)) + B)
    xlim([-xinit,xinit])
    ylim([-4,4])
    hold off
    pause(plotter.pause)
    
    if D<0
        hold on; title("collision"); hold off
    end
end
end