clear all
close all
n=9;
A_star=optimvar('As',3,3);
A_zero=optimvar('Az',3,n);
E_star_zero = optimvar('E',n,3);
B=optimvar('B',3,2);
load('./../DataSet/DataSet_train.mat');
fun=(A_star*Xt_1)+(A_zero*tanh(E_star_zero*Xt_1))+(B*Ut_1);
obj = sum(sum(transpose((fun - Xt).^2)));
lsqproblem = optimproblem("Objective",obj);
X0.As=rand(3,3)*0.01;
X0.Az=rand(3,n)*0.01;
X0.B=eye(3,2)./154;
X0.E=rand(n,3)*0.01;
options=optimoptions(lsqproblem,'Algorithm','levenberg-marquardt','Display','iter','MaxIterations',400)
show(lsqproblem)
[sol,fval] = solve(lsqproblem,X0,'Options',options)

O_As=sol.As;
O_Az=sol.Az;
O_B=sol.B;
O_E=sol.E;

save('Optimom_final_tanh_9','O_Az','O_As','O_B','O_E')