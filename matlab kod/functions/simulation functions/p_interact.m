function pea = p_interact(S,beta_EA,RUprop)
% computes the probability that driver will react based on regression
% parameters beta=[bet0,beta1,beta2].
A = S(1).pos;
B = S(2).pos;
speedA = S(1).speed;
speedB = S(2).speed;
beta = [beta_EA.int;beta_EA.dist;beta_EA.speed];

D = norm(A-B);
X1 = [1, D, speedB];
X2 = [1, D, speedA];
X = [X1;X2];
pea = RUprop.alert./(1 + exp(X*beta));
end