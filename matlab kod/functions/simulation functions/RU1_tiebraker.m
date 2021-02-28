function RU1 = RU1_tiebraker(state,r)
% function that calculates which road user is enters conflict region first
% in case there is a tie

% extract states
stateA = state(1);
stateB = state(2);

% get parameters for lines that define boundaries of tangent path
lineparA = find_path_boundaries(stateA,r);
lineparB = find_path_boundaries(stateB,r);

% compute time t when they first touch boundaries of others path
t1A = time_to_line(stateA, lineparB(1,:),r);
t2A = time_to_line(stateA, lineparB(2,:),r);
tA = min(t1A,t2A);
t1B = time_to_line(stateB, lineparA(1,:),r);
t2B = time_to_line(stateB, lineparA(2,:),r);
tB  = min(t1B,t2B);

[~, RU1] = min([tA,tB]);


% A0 = stateA.pos - stateA.speed*exp(1i*stateA.theta);
% A1 = stateA.pos;

% B0 = stateB.pos - stateB.speed*exp(1i*stateB.theta);
% B1 = stateB.pos;
% hold on 
% plot_circle(A1,r,"black")
% hold on
% plot_circle(B1,r,'black')

% plot_circle(A0 + tA*(A1-A0),r,'green')
% plot_circle(B0 + tB*(B1-B0),r,[255,17,227]/255)
% 
% drawline(lineparA(1,:),[43,163,255]/255)
% drawline(lineparA(2,:),'blue')
% drawline(lineparB(1,:),[255,106,0]/255)
% drawline(lineparB(2,:),'red')
% plot_circle(A0,r,'cyan')
% plot_circle(B0,r,'cyan')


%%
end
