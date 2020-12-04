function plot_pos(S, plotter, xinit, r, int_state) 
% plots position of each disc based on their positions in state S.

% unpack variables
A = S(1).pos;
B = S(2).pos;
reactionA = int_state.a;
reactionB = int_state.b;

if plotter.enc == 1
    xlim([-xinit,xinit])
    ylim([-4,4])
    
    plot((r*cos(linspace(0,2*pi,50))+real(A))+1i*(r*sin(linspace(0,2*pi,50)) + imag(A)))
    title(sprintf("collision, reactionA = %d, reactionB = %d", reactionA,reactionB))
    hold on
    
    plot((r*cos(linspace(0,2*pi,50))+real(B))+1i*(r*sin(linspace(0,2*pi,50))+imag(B)))
    xlim([-xinit,xinit])
    ylim([-4,4])
    hold off
    pause(plotter.pause)
end
end