clear all
close all


p=50;


Num_of_Signals=1;
X0=[0;0;0];
out=sim('Model.slx');
u=out.input.data(:,1,:);
u=reshape(u,[2 6001]);
u=u(:,1:6000);
x=out.output.data(:,1,:);
x=reshape(x,[3 6001]);
xt_1=x(:,1:6000);
xt=x(:,2:6001);

Ut_1=u;
Xt_1=xt_1;
Xt=xt;

f = waitbar(0,'Generating...');
for h = 1:3
    Num_of_Signals=h;
    for c = 1:p
        waitbar(((h-1)*p+c)/((3*p)+15),f,sprintf('Generating your data %d of %d',((h-1)*p+c),3*p+15));
        
        X0=[rand*100;rand*100;rand*100];
        out=sim('Model.slx');

        u=out.input.data(:,1,:);
        u=reshape(u,[2 6001]);
        u=u(:,1:6000);
        x=out.output.data(:,1,:);
        x=reshape(x,[3 6001]);
        xt_1=x(:,1:6000);
        xt=x(:,2:6001);

        Ut_1=cat(2,Ut_1,u);
        Xt_1=cat(2,Xt_1,xt_1);
        Xt=cat(2,Xt,xt);



    end
end
X0=[0;0;0];
for h = 1:3
    Num_of_Signals=h;
    for c = 1:5
        waitbar(((h-1)*5+c)/(3*p+15),f,sprintf('Generating your data %d of %d',+((h-1)*5+c),(3*p+15)));

        out=sim('Model.slx');

        u=out.input.data(:,1,:);
        u=reshape(u,[2 6001]);
        u=u(:,1:6000);
        x=out.output.data(:,1,:);
        x=reshape(x,[3 6001]);
        xt_1=x(:,1:6000);
        xt=x(:,2:6001);

        Ut_1=cat(2,Ut_1,u);
        Xt_1=cat(2,Xt_1,xt_1);
        Xt=cat(2,Xt,xt);



    end
end
close(f)

number_of_samples=size(Ut_1);
number_of_samples=number_of_samples(2);
s = randperm(number_of_samples);
Ut_1(1,:)=Ut_1(1,s);
Ut_1(2,:)=Ut_1(2,s);

Xt_1(1,:)=Xt_1(1,s);
Xt_1(2,:)=Xt_1(2,s);
Xt_1(3,:)=Xt_1(3,s);

Xt(1,:)=Xt(1,s);
Xt(2,:)=Xt(2,s);
Xt(3,:)=Xt(3,s);


disp('Size of U(t-1):')
disp(size(Ut_1))
disp('Size of X(t-1):')
disp(size(Xt_1))
disp('Size of X(t):')
disp(size(Xt))

save('DataSet_train','Ut_1','Xt_1','Xt')
disp('DataSet has been saved to DataSet_train.mat file')