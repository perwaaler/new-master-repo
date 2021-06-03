function out = path_analysis(ind,stateA,stateB,opt,r,plotting)
% k is the total path length. k is number of steps taken, or how far into
% the future the path is generated.
k         = 80;
% path2.pos = path2.pos(ind);
path{1}  = gener_pred_path_iter(stateA, k,ind);
path{2}  = gener_pred_path_iter(stateB, k,ind);
k1     = path{1}.length;
k2     = path{2}.length;
tstop1 = path{1}.tstop;
tstop2 = path{2}.tstop;

dists = abs(path{1}.padded_pos - path{2}.padded_pos);
loc_extr = find_vector_min(dists);

%%% translate complex coordinates to cartesian coordinates
x1 = real(path{1}.pos);
y1 = imag(path{1}.pos);

x2 = real(path{2}.pos);
y2 = imag(path{2}.pos);

% parameter values (time)
t1 = path{1}.ind - 1;
t2 = path{2}.ind - 1;

% corresponding output values (2d-positions)
v1 = [x1;y1];
v2 = [x2;y2];

% special case when the RU is standing still
if size(v1,2)==1
    v1 = [v1,v1]; %#ok<*AGROW>
    t1 = [0,1];
end
if size(v2,2)==1
    v2 = [v2,v2];
    t2 = [0,1];
end


% fit splines
splineA = makima(t1,v1);
splineB = makima(t2,v2);
ppA =@(t) ppval(splineA,t);
ppB =@(t) ppval(splineB,t);
% create piecewise fcn which ensures no crazy behviour for RUs that have
% stopped
ppA =@(t) pw_fcn(ppA,t,tstop1);
ppB =@(t) pw_fcn(ppB,t,tstop2);

% prevent algorithm from selecting times outside of the generated path if
% the acceleration is negative
tstop1 = min(k1,tstop1);
tstop2 = min(k2,tstop2);

%%%% functions that determines search range of fminsearch:
% discourage optimum after RU stops
bound_stop = @(t) (1+exp((t(1)-tstop1)/0.001))*...
                  (1+exp((t(2)-tstop2)/0.001));
% fcn that discourages fmin to search for negative minima
bound_noneg = @(t) (1+exp(-t/0.001));

% define distance functions
d_prox_t =@(t) norm(ppA(t(1)) - ppB(t(2)))*...
                bound_noneg(t(1))*bound_noneg(t(2))*...
                bound_stop(t);
d_min_t  =@(t) norm(ppA(t) - ppB(t))*...
                bound_noneg(t);

%%%% optimization %%%%
% find proximity points, minimum sep. points, and corresponding times

[lbnd,rbnd] = search_range_dmin(loc_extr, ind);
[t_min,d_min] = fminbnd(d_min_t, lbnd, rbnd, opt);
% opt = optimset('TolX',0.01,'MaxIter',30,'Display','off');
% [t_min(1),d_min(1)] = fminbnd(d_min_t, 0, k/2, opt);
% [t_min(2),d_min(2)] = fminbnd(d_min_t, k/2,k, opt);
% [t_min,I] = min(t_min);
% d_min = d_min(I);

[t_prox,d_prox] = fminsearch(d_prox_t, [t_min,t_min], opt);
% [t_prox(1,:),d_prox(1)] = fminsearch(d_prox_t, [t_min,t_min], opt);
% t_prox = t_prox(1,:);
% [t_prox(2,:),d_prox(2)] = fminsearch(d_prox_t, flip(t_prox(1,:)), opt);
% [d_prox, I_global] = min(d_prox);
% t_prox = t_prox(I_global,:);

if d_prox < 0.6
    [t_prox(1), d_prox(1)] = findTTO(d_prox_t, t_prox, r, 1,[k1,k2],opt);
    [t_prox(2), d_prox(2)] = findTTO(d_prox_t, t_prox, r, 2,[k1,k2],opt);
end

% find tangent lines at minimum points
h = 0.0001;
slope_minA  = find_slope(path{1},t_min,ppA,h); 
slope_minB  = find_slope(path{2},t_min,ppB,h); 

% compute proximity and minimum points
Aprox = cart2complex(ppA(t_prox(1)));
Bprox = cart2complex(ppB(t_prox(2)));
Amin  = cart2complex(ppA(t_min));
Bmin  = cart2complex(ppB(t_min));

Amin_pert  = Amin + 1 + 1i*slope_minA;
Bmin_pert  = Bmin + 1 + 1i*slope_minB;

orient(1) = find_side(1 + 1i*slope_minA,Bmin-Amin);
orient(2) = find_side(1 + 1i*slope_minB,Amin-Bmin);

out.orient = orient;
out.t_prox = t_prox;
out.t_min  = t_min;
out.d_prox = d_prox;
out.d_min  = d_min;

%%%% plotting %%%%
% points to evaluate function (for plotting only)
if plotting==1
    
    % plot generated path
    hold on
    for i=1:length(ind)
        a = plot_circle(path{1}.padded_pos(i),0.3,"green");
        if i==length(path{1}.ind)
            a.Color = "yellow";
        end
        a = plot_circle(path{2}.padded_pos(i),0.3,"cyan");
        if i==length(ind)
            a.Color = "yellow";
        end
    end
%     color starting points red
    a = plot_circle(path{1}.pos(1), 0.3, col2trip("grey"));
    a.LineWidth = 1;
    b = plot_circle(path{2}.pos(1), 0.3, col2trip("grey"));
    b.LineWidth = 1;

    
    tt1 = linspace(0,k1,200);
    tt2 = linspace(0,k2,200);
    V1 = ppval(splineA, tt1);
    V2 = ppval(splineB, tt2);

    % plot spline curve fits
    a=plot(V1(1,:),V1(2,:),'-','color','green');
    a.LineWidth=1;
    a=plot(V2(1,:),V2(2,:),'-','color','blue');
    a.LineWidth=1;

    % plot proximity points
    a=plot_circle(Amin,.3,"black");
    a.LineWidth=1;
    b=plot_circle(Bmin,.3,"black");
    b.LineWidth=1;
    drawline([Amin, Amin_pert],'black');
    drawline([Bmin, Bmin_pert],'black');
    title(sprintf("(%s,%s), t_{min} = %0.01f, t_{prox}=(%0.01f,%0.01f)",...
                    orient(1),orient(2),t_min,t_prox(1),t_prox(2)))

    xlim([-10,10])
    a = plot_circle(Aprox,.3,col2trip("darkgreen"));
    a.LineWidth = 2;
    b = plot_circle(Bprox,.3,col2trip("darkblue"));
    b.LineWidth = 2;

    pause(1)
end

end