clear all
close all

load('Optimom_final_tanh_11.mat')
load('./../DataSet/DataSet_test.mat')
%0...1...2....3
a=size(Ut_1);
for j=0:((a(2)/6000)-1)
    t=j*6000;
    xt_1=Xt_1(:,t+1:t+1);
    Xt_hat=zeros(3,6000);
    for i=1:6000

        Xt_hat(:,i)=(O_As*xt_1)+(O_Az*tanh(O_E*xt_1))+(O_B*Ut_1(:,i+t));
        xt_1=Xt_hat(:,i);

    end
    figure
    p=plot(transpose(Xt(:,t+1:t+6000)));
    p(1).DisplayName="x1";
    p(2).DisplayName="x2";
    p(3).DisplayName="x3";
    hold on
    p=plot(transpose(Xt_hat));
    p(1).LineStyle="--";
    p(1).DisplayName="x1 Estimation";
    p(2).LineStyle="--";
    p(2).DisplayName="x2 Estimation";
    p(3).LineStyle="--";
    p(3).DisplayName="x3 Estimation";
    legend
    title(sprintf('System States and their Estimations (Test %d)',j+1))
    xlabel('Sample')
    
    figure
    p=plot(transpose(Xt_hat)-transpose(Xt(:,t+1:t+6000)));
    p(1).DisplayName="x1 Error";
    p(2).DisplayName="x2 Error";
    p(3).DisplayName="x3 Error";
    title(sprintf('System States and their Estimations Error (Test %d)',j+1))
    xlabel('Sample')
    legend
end