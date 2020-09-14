function plot_pos(A0, B0, pause_length, xinit, r, plot_setting,reactionA,reactionB)            
if plot_setting == 1
    xlim([-xinit,xinit])
    ylim([-4,4])
    plot([r*cos(linspace(0,2*pi,50))+real(A0)]+1i*[r*sin(linspace(0,2*pi,50))+imag(A0)])
    title(sprintf("collision, reactionA = %d, reactionB = %d", reactionA,reactionB))
    hold on
    plot([r*cos(linspace(0,2*pi,50))+real(B0)]+1i*[r*sin(linspace(0,2*pi,50))+imag(B0)])
    xlim([-xinit,xinit])
    ylim([-4,4])
    hold off
    pause(pause_length)
end
end