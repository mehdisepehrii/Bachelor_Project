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
Ka=inv(value(U))*value(N_bar);
out=sim("Fault_Estimator.slx");
x_hat=zeros(4,6001);
Fa_hat=zeros(1,6001);
Y=reshape(out.y.data,[3,6001]);
X=reshape(out.X.data,[3,6001]);
Fs=reshape(out.Fs.data,[1,6001]);
Fa=reshape(out.Fa.data,[1,6001]);
u=reshape(out.u.data,[2,6001]);
x_hat(:,1)=[Y(:,1)' 0]';
for i=2:6001
    if mod(i,10)==0
        theta=(tanh(O_E*x_hat(1:3,i-1)))./(O_E*x_hat(1:3,i-1));
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
        Ka=inv(value(U))*value(N_bar);
    end
    x_hat(:,i)=(A_bar*x_hat(:,i-1))+(B_bar*u(:,i-1))+(G*A_zero*(tanh(E_zero*x_hat(:,i-1))))+(L_bar*Y(:,i))+(Ka*(Y(:,i-1)-(C*x_hat(:,i-1))));
    Fa_hat(:,i-1)=(2/3)*H*((eye(p)-C*N)*Y(:,i)-C*T*A*x_hat(:,i-1)-C*T*B*u(:,i-1)-C*T*A_zero*(tanh(E_zero*x_hat(:,i-1))));
end

plot([X',transpose(x_hat(1:3,:))])
figure
plot([Fs',transpose(x_hat(4,:))])
figure
plot([Fa' Fa_hat'])