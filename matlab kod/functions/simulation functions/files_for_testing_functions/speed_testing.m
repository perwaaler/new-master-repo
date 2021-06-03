% testing different ways to generate path

j=20;
tic
pos1 = gener_pred_path(state,j);
toc
tic
gener_pred_path_iter(state,j);
toc

% clf
% plot_circle(state.pos,.3)
% for i=1:k
%     plot_circle(pos(i),.3)
% end

tic
A0 = state.pos;
speed0 = state.speed;
dspeed = state.dspeed;
theta0 = state.theta;
dtheta = state.dtheta;

pos = zeros(1,j);
pos(1) = A0;
for i=1:j-1
    pos(i+1) = pos(i) + (state.speed + i*state.dspeed)*...
        exp(1i*(state.theta + i*state.dtheta));
end
pos2 = [A0,pos];
toc



% for k = 80, the function that uses matrix multiplication is 
% twice as fast. For large k, the for-loop becomes faster.


%% effect of combining for loops into one
   



n = 10^3;
x = cell(1,n);
y = cell(1,n);
n_tsample = 400;
time_frac = zeros(1,n_tsample);
t1 = zeros(1,n_tsample);
t2 = zeros(1,n_tsample);
for j=1:n_tsample
    tic
    for i=1:n
        X = rand(10,10);
        x{i} = X*X;
        Y = rand(10,10);
        y{i} = Y*Y;
    end
    
    tic
    for i=1:n
        X = rand(10,10);
        x{i} = X*X;
    end
    for i=1:n
        Y = rand(10,10);
        y{i} = Y*Y;
    end
    
    % one for loop

    t1(j) = toc;

    t2(j) = toc; 
end




time_frac = t1./t2;


mean(time_frac) + 1.96*std(time_frac)*[-1,0,1]/sqrt(n_tsample) %#ok<NOPTS>
clf
plot(time_frac)
% There seems to be little to no difference in performance, judging by the
% confidence intervals. Be careful! If the for loops are run separately
% this can dramatically influence the test-results. Run them together, and
% try changing order to be certain.
%%



tic
X = ones(10*n,10);
for i=1:n
    y{i} = X((i-1)*10 + 1 :i*10,:)*X((i-1)*10 + 1 :i*10,:);
end
toc
tic
for i=1:n
    X = ones(10,10);
    x{i} = X*X;
end
toc

% it appears that constructing X inside the forloop instead of creating a
% large matrix before the for loop is far more efficient. The latter loop
% is more than twice as fast.
%%
state1.pos = normrnd(-2,1) + normrnd(-2.5,3)*1i;
state1.theta = normrnd(deg2rad(10),deg2rad(7));
state1.dtheta = normrnd(deg2rad(10),deg2rad(6));
state1.speed = gam_rnd(0.16,0.02);
state1.dspeed = normrnd(0.01,0.01);  
state1.RUprop = RUprop;
state1.id     = 1;

n=300;
time_ratio = zeros(1,n);
asave = zeros(1,n);
bsave = zeros(1,n);
for i=1:n
tic
k = 25; % number of steps taken. Total path lenght is k+1
gener_pred_path_iter(state1,max_delta(1),k);
a = toc;
asave(i)=a;
tic
k = 100; % number of steps taken. Total path lenght is k+1
gener_pred_path_iter(state1,max_delta(1),k);
b = toc;
bsave(i)=b;

time_ratio(i) = a/b;
end

clf
plot(asave)


for k=1:100
g =@() gener_pred_path_iter(state1,max_delta(1),500);
bsave(k) = timeit(g);
f =@() gener_pred_path_iter(state1,max_delta(1),100);
asave(k) = timeit(f);
% end
% for k=1:100
end

mean(bsave)/mean(asave)


[path1,slopes1,K1] = gener_pred_path_iter(state1,max_delta(1),k);



%% testing the new step algorithm
% states = gen_random_states();
% stateA = states{1};
% stateB = states{2};




