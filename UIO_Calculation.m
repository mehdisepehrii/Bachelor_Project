clear all
close all

load('./Learning/Optimom_final_tanh_7.mat')

n=3; %number of system states
m=2; %number of system inputs
p=3; %number of system output
r=7; %number of neurons 
s=1; %number of sensors faults
q=1; %number of actuators faults 

C_star=eye(p,n);
F=[0;0;1]; %Fault Can only apear on Sensor x2
W=zeros(n,n);%Disturbance Matrix
E=[eye(n) zeros(n,s)];
E_zero=O_E*E;
A=[O_As zeros(n,s)];
C=[C_star F];
temp=[1,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0;0,0,-1,0,0,1];
T=temp(1:n+s,1:n);
N=temp(1:n+s,n+1:n+p);
B=O_B;
D=[1/154;0;0];%Actuator Fault Matrix
A_zero=O_Az;
H=inv(transpose(C*T*D)*C*T*D)*transpose(C*T*D);
G=T-T*D*H*C*T;
W_bar=G*W;
A_bar=G*A;
B_bar=G*B;
L_bar=(N-T*D*H*C*N+T*D*H);

theta=(tanh(O_E*[20;60;15]))./(O_E*[20;60;15]);
temp=A_zero*diag(theta)*E_zero;
A_tilda_alpha=G*(A+temp);
A_3_alpha=A+temp;
H1=T'*C'*H'*H*C*T;
mu=sdpvar(1);
P=sdpvar(n+s,n+s,'full');
U=sdpvar(n+s,n+s,'full');
N_bar=sdpvar(n+s,p,'full');
U_A_2_alpha=U*A_tilda_alpha-N_bar*C;
lmi=[A_3_alpha'*H1*A_3_alpha-P,A_3_alpha'*H1*W,U_A_2_alpha';W'*H1*A_3_alpha,W'*H1*W-mu^2*eye(n),W_bar'*U';U_A_2_alpha,U*W_bar,P-U-U'];
F = [lmi<=0; P>=0 ; mu>=0;N_bar<=1];
optimize(F, mu);
mu=value(mu)
Ka=inv(value(U))*value(N_bar)

save('Fault_Estimator_Parameters.mat','G','C','p','L_bar','A_bar','B_bar','A_zero','E_zero','C','T','B','A','N','H','Ka')