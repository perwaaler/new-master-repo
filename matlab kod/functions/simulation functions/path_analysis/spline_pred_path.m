function Tsep = spline_pred_path(stateA,stateB,k,sparse_ind,genInfo)
% k is the total path length. k is number of steps taken, or how far into
% the future the path is generated.

if nargin==4
    genInfo.nIter = 5;
    genInfo.genPath = false;
    genInfo.analyze = false;
    genInfo.plotRes = false;
end

r = stateA.RUprop.r;
paths = gener_pred_path_iter_test(stateA,stateB,k,sparse_ind,genInfo);
pathA = paths.A;
pathB = paths.B;

if genInfo.plotRes
    clf
    hold on
    for i=1:length(pathA.full)
        plot(pathA.full(i),'Marker','.','Color','green');
        if ismember(i,pathA.s_ind)
            plot_circle(pathA.full(i),.3,col2trip('lightgreen'));
        end
    end
    for i=1:length(pathB.full)
        plot(pathB.full(i),'Marker','.','Color','blue');
        if ismember(i,pathB.s_ind)
            plot_circle(pathB.full(i),.3,col2trip('lightblue'));
        end
    end
    % color starting points red
    plot_circle(pathA.pos_sparse(1), 0.3, col2trip("grey"),4);
    plot_circle(pathB.pos_sparse(1), 0.3, col2trip("grey"),4);
    plot_circle(pathA.full(paths.t_prox(1)+1), 0.3, col2trip('lightgreen'),4);
    plot_circle(pathB.full(paths.t_prox(2)+1), 0.3, col2trip('lightblue'),4);
    plot_circle(pathA.full(paths.t_min+1), 0.3, col2trip("lightgrey"),4);
    plot_circle(pathB.full(paths.t_min+1), 0.3, col2trip("lightgrey"),4);
end

%%% translate complex coordinates to cartesian coordinates
x1 = real(pathA.pos_sparse);
y1 = imag(pathA.pos_sparse);

x2 = real(pathB.pos_sparse);
y2 = imag(pathB.pos_sparse);

% parameter values (time)
t1 = ind2time(pathA.s_ind);
t2 = ind2time(pathB.s_ind);

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

t_relA = paths.rel_time(1);
t_relB = paths.rel_time(2);
ppA =@(t) pw_fcn(ppA,t,t_relA);
ppB =@(t) pw_fcn(ppB,t,t_relB);

% prevent algorithm from selecting times outside of the generated path if
% the acceleration is negative

%%%% functions that determines search range of fminsearch:
% discourage optimum after RU stops
bound_stop = @(t) (1+0.2*exp((t(1)-t_relA)/0.01))*...
                  (1+0.2*exp((t(2)-t_relB)/0.01));
% fcn that discourages fmin to search for negative minima
bound_noneg = @(t) (1+0.2*exp(-t/0.01));

% define distance functions
d_prox_t =@(t) norm(ppA(t(1)) - ppB(t(2)))*...
                bound_noneg(t(1))*bound_noneg(t(2))*...
                bound_stop(t);
d_min_t  =@(t) abs(norm(ppA(t) - ppB(t))-2*r)*...
                bound_noneg(t);

            

%%%% optimization %%%%
% find proximity points, minimum sep. points, and corresponding times
opt = optimset('TolX',0.01,'MaxIter',genInfo.nIter,'Display','off');

tmin0 = paths.t_min;
lbnd = max(0,tmin0 - 1);
ubnd = min(  tmin0 + 1,max(t_relA,t_relB));
[t_min,d_min] = fminbnd(d_min_t, lbnd, ubnd, opt);

t_prox=[0,0];
if paths.contact
    [t_prox(1), dproxA] = findTTO_V2(d_prox_t, paths, r, 1,opt);
    [t_prox(2), dproxB] = findTTO_V2(d_prox_t, paths, r, 2,opt);
else
    [t_prox, dproxA] = findTTO(d_prox_t, paths, r, 1,opt);
    dproxB = dproxA;
end    

% find tangent lines at proximity and minimum points
h = 0.0001;
angle_minA  = find_angle(paths,1,t_min,ppA,h); 
angle_minB  = find_angle(paths,2,t_min,ppB,h); 

% calculate minimality points and differentials to find orientation
Amin  = cart2complex(ppA(t_min));
Bmin  = cart2complex(ppB(t_min));
dAmin = exp(1i*angle_minA);
dBmin = exp(1i*angle_minB);

orient(1) = find_side(dAmin, Bmin-Amin);
orient(2) = find_side(dBmin, Amin-Bmin);
phi = find_action_angle(orient);

Tsep.t_prox = t_prox;
Tsep.t_min  = t_min;
Tsep.d_prox = [dproxA,dproxB];
Tsep.d_min  = d_min;
Tsep.actionAngle = phi;



%%% plotting %%%%
if genInfo.plotRes
% hold on
% for i=1:length(pathA.full)
%     a = plot_circle(pathA.full(i),0.3,"green");
%     if i==length(pathA.s_ind)
%         a.Color = "yellow";
%         a.LineWidth = 2;
%     end
% end
% for i=1:length(pathB.full)
%     a = plot_circle(pathB.full(i),0.3,"cyan");
%     if i==length(pathB.s_ind)
%         a.Color = "yellow";
%         a.LineWidth = 2;
%     end
% end


tt1 = linspace(0,t_relA,200);
tt2 = linspace(0,t_relB,200);
V1 = ppval(splineA, tt1);
V2 = ppval(splineB, tt2);

% plot spline curve fits
plot(V1(1,:),V1(2,:),'-','color','green','LineWidth',1);
plot(V2(1,:),V2(2,:),'-','color','blue' ,'LineWidth',1);

% compute proximity points and minimality points
Aprox = cart2complex(ppA(t_prox(1)));
Bprox = cart2complex(ppB(t_prox(2)));
% comopute how to steer to maximize d_min
action = orient2action(orient);

% draw minimum positions
plot_circle(Amin,.3,"black",2);
plot_circle(Bmin,.3,"black",2);
drawline([Amin, Amin + dAmin],'black');
drawline([Bmin, Bmin + dBmin],'black');
title(sprintf("A %s,B %s, t_{min} = %0.01f, t_{prox}=(%0.01f,%0.01f)",...
                            action(1),action(2),t_min,t_prox(1),t_prox(2)))

xlim([-10,10])

% draw proximity points
plot_circle(Aprox,.3,col2trip("darkgreen"),2);
plot_circle(Bprox,.3,col2trip("darkblue"),2);
pause(genInfo.plotP)
end

end
