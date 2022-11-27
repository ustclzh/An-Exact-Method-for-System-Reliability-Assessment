%作出G=\alpha 的图像（三维曲面）



clear;
clc;
n_1=10;
n_2=10;
n_3=10;
global N;
N=10000;

%%
%参数设定
%所有随机数在这里生成
global randnum_1 randnum_2 randnum_31 randnum_32 Q_1 Q_2 Q_31 Q_32
randnum_1=normrnd(0,1,2*n_1,N);
randnum_2=rand(n_2,N);
randnum_31=normrnd(0,1,1,N);
randnum_32=normrnd(0,1,n_3-1,N);
Q_1=normrnd(0,1,n_1,N);
Q_2=rand(n_2,N);
Q_31=normrnd(0,1,1,N);
Q_32=normrnd(0,1,n_3-1,N);
aa=0.2; %置信水平alpha
eta=150;
m=1.5;
theta=240;
sigma=0.8;
mu=4;
t_0=20;
c_real=exp(-t_0/theta)*exp(-(t_0/eta)^m)*normcdf((mu-log(t_0))/sigma);%真值
C_real=[exp(-t_0/theta);exp(-(t_0/eta)^m);normcdf((mu-log(t_0))/sigma)];


%%
%在三维空间中搜索G=\alpha的点
    %生成一组数据
    data_1=-theta*log(rand(1,n_1));%指数
    data_2=eta*(-log(rand(1,n_2))).^(1/m);%威布尔
    data_3=exp(normrnd(mu,sigma,1,n_3));%对数正态
    data=[data_1;data_2;data_3];
    data_2=log(data_2);
    data_3=log(data_3);
    r_hat_1=exp(-t_0/(mean(data_1)));
    r_hat_2=exp(-exp((log(t_0)-mean(data_2))/(6^0.5*((n_2-1)*var(data_2)/n_2)^0.5/pi)-psi(1)));
    r_hat_3=normcdf((mean(data_3)-log(t_0))/(((n_3-1)*var(data_3)/n_3)^0.5));
    c=r_hat_1*r_hat_2*r_hat_3;%系统可靠性的估计量

    a_u=1;
    a_l=0;
    while a_u-a_l>0.0000001
    a1=(a_u+a_l)/2;
        if G(a1,1,1,n_1,n_2,n_3,c,aa)<=0
            a_l=a1;
            a1=0.5*(a1+a_u);
        end
        if G(a1,1,1,n_1,n_2,n_3,c,aa)>=0
            a_u=a1;
            a1=0.5*(a1+a_l);
        end
    end
    a_u=1;
    a_l=0;
    while a_u-a_l>0.0000001
    a2=(a_u+a_l)/2;
        if G(1,a2,1,n_1,n_2,n_3,c,aa)<=0
            a_l=a2;
            a2=0.5*(a2+a_u);
        end
        if G(1,a2,1,n_1,n_2,n_3,c,aa)>=0
            a_u=a2;
            a2=0.5*(a2+a_l);
        end
    end
    a_u=1;
    a_l=0;
    while a_u-a_l>0.0000001
    a3=(a_u+a_l)/2;
        if G(1,1,a3,n_1,n_2,n_3,c,aa)<=0
            a_l=a3;
            a3=0.5*(a3+a_u);
        end
        if G(1,1,a3,n_1,n_2,n_3,c,aa)>=0
            a_u=a3;
            a3=0.5*(a3+a_l);
        end
    end
a=[a1,a2,a3];
A=zeros(10000,3);
for i=1:100
    i
    for j=1:100
        u1=a1+i*(1-a1)/100;
        u2=a2+i*(1-a2)/100;
        r_u=1;
        r_l=a3;
        while r_u-r_l>0.0001;
            r=(r_u+r_l)/2;
            if G(u1,u2,r,n_1,n_2,n_3,c,aa)<=0
            r_l=r;
            r=0.5*(r+r_u);
            end
            if G(u1,u2,r,n_1,n_2,n_3,c,aa)>=0
            r_u=r;
            r=0.5*(r+r_l);
            end
        end
A(100*(i-1)+j,:)=[u1,u2,r];
    end
end

    
    
    