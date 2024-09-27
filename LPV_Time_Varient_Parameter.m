clear all
close all
n=10;
addpath('Learning')
load('./DataSet/DataSet_train.mat');
load('./Learning/Optimom_final_2.mat')
A_star=O_As;
A_zero=O_Az;
B=O_B;
E_star_zero=O_E;

theta=optimvar('theta',n)
fun2=(A_star*Xt_1)+(A_zero*diag(theta)*(E_star_zero*Xt_1))+(B*Ut_1);
obj = sum(sum(transpose((fun2 - Xt).^2)));
lsqproblem = optimproblem("Objective",obj);
X0.theta=ones(n,1)./10;
options=optimoptions(lsqproblem,'Display','iter','MaxIterations',400)
show(lsqproblem)
[sol,fval] = solve(lsqproblem,X0,'Options',options)

theta=sol.theta;

save('theta','theta')