
clear;
clc;
%%
n_1=2;
n_2=3;
n_3=4;
global N;
N=10000;
%所有随机数在这里生成
global randnum_1 randnum_2 randnum_3 Q_1 Q_2 Q_3 
randnum_1=normrnd(0,1,2*n_1,N);%用于生成自由度为2n_1的卡方分布
randnum_2=normrnd(0,1,2*n_2,N);%用于生成自由度为2n_1的卡方分布
randnum_3=normrnd(0,1,2*n_3,N);%用于生成自由度为2n_1的卡方分布

Q_1=normrnd(0,1,n_1,N);
Q_2=normrnd(0,1,n_2,N);
Q_3=normrnd(0,1,n_3,N);

aa=0.1; %置信水平alpha
%%    
%数据输入，直接给出估计量
r_hat_1=exp(-1/17.9855);
r_hat_2=exp(-1/7.305);
r_hat_3=exp(-1/31.271);
c=r_hat_1*r_hat_2*r_hat_3;%系统可靠性的估计量
%求一个初始值
r=0.5;
r_u=1;
r_l=0;
    while  sum(r_u-r_l)>0.000001
        if G(1,1,r,n_1,n_2,n_3,c,aa)<0
            r_l=r;
            r=0.5*(r+r_u);
        end
        if G(1,1,r,n_1,n_2,n_3,c,aa)>0
            r_u=r;
            r=0.5*(r+r_l);
        end
        if G(1,1,r,n_1,n_2,n_3,c,aa)==0
            break
        end

    end
%得到一个可行解:
r=[1,1,r];


%%
%算法参数调整
beta=0.05;
lambda=1;
alpha=0.4;%初始步长
%算法参数调整
%Lagrange Newton算法
P=1;
P_temp=0;
%while abs(P-P_temp)>0.00001
for iiii=1:30
    
    P_temp=P;
    P=(r(2)*r(3)-lambda*f_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N))^2+(r(1)*r(3)-lambda*f_2(r(1),r(2),r(3),n_1,n_2,n_3,c,N))^2+(r(1)*r(2)-lambda*f_3(r(1),r(2),r(3),n_1,n_2,n_3,c,N))^2+G(r(1),r(2),r(3),n_1,n_2,n_3,c,aa)^2;
    
    f12=(f_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N)-f_1(r(1),r(2)-0.01,r(3),n_1,n_2,n_3,c,N))/0.01;%二阶偏导数，用差商代替
    f13=(f_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N)-f_1(r(1),r(2),r(3)-0.01,n_1,n_2,n_3,c,N))/0.01;
    f23=(f_2(r(1),r(2),r(3),n_1,n_2,n_3,c,N)-f_2(r(1),r(2),r(3)-0.01,n_1,n_2,n_3,c,N))/0.01;
    
    W=[-lambda*df_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N),r(3)-lambda*f12,r(2)-lambda*f13,-f_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N);
        r(3)-lambda*f12,-lambda*df_2(r(1),r(2),r(3),n_1,n_2,n_3,c,N),r(1)-lambda*f23,-f_2(r(1),r(2),r(3),n_1,n_2,n_3,c,N);
        r(2)-lambda*f13,r(1)-lambda*f23,-lambda*df_3(r(1),r(2),r(3),n_1,n_2,n_3,c,N),-f_3(r(1),r(2),r(3),n_1,n_2,n_3,c,N);
        -f_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N),-f_2(r(1),r(2),r(3),n_1,n_2,n_3,c,N),-f_3(r(1),r(2),r(3),n_1,n_2,n_3,c,N),0];         
    
    y=-[r(2)*r(3)-lambda*f_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N);r(1)*r(3)-lambda*f_2(r(1),r(2),r(3),n_1,n_2,n_3,c,N);r(1)*r(2)-lambda*f_3(r(1),r(2),r(3),n_1,n_2,n_3,c,N);-G(r(1),r(2),r(3),n_1,n_2,n_3,c,aa)];
    
    x=W\y;
    lambda_1=lambda+alpha*x(4);
    r_temp=r+alpha*x(1:3)';
    if max(r_temp)>1
        alpha=alpha/2;
        continue;
    end
    if min(r_temp)<0
        alpha=alpha/2;
        continue;
    end
    P_1=(r_temp(2)*r_temp(3)-lambda_1*f_1(r_temp(1),r_temp(2),r_temp(3),n_1,n_2,n_3,c,N))^2+(r_temp(1)*r_temp(3)-lambda_1*f_2(r_temp(1),r_temp(2),r_temp(3),n_1,n_2,n_3,c,N))^2+(r_temp(1)*r_temp(2)-lambda_1*f_3(r_temp(1),r_temp(2),r_temp(3),n_1,n_2,n_3,c,N))^2+G(r_temp(1),r_temp(2),r_temp(3),n_1,n_2,n_3,c,aa)^2;
    if P_1<=(1-beta*alpha)*P
        r=r+alpha*x(1:3)';
        lambda=lambda+alpha*x(4);
    else
        alpha=alpha/2;
    end
  R_L=r(1)*r(2)*r(3);
end
%Lagrange Newton算法
%%




