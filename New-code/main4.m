clear;
clc;
n_1=10;
n_2=10;
n_3=10;
global N;
N=10000;


%�������������������
global randnum_1 randnum_2 randnum_31 randnum_32 Q_1 Q_2 Q_31 Q_32
randnum_1=normrnd(0,1,2*n_1,N);
randnum_2=rand(n_2,N);
randnum_31=normrnd(0,1,1,N);
randnum_32=normrnd(0,1,n_3-1,N);
Q_1=normrnd(0,1,n_1,N);
Q_2=rand(n_2,N);
Q_31=normrnd(0,1,1,N);
Q_32=normrnd(0,1,n_3-1,N);


aa=0.1; %����ˮƽalpha
eta=180;
m=1.5;
theta=220;
sigma=0.5;
mu=4;
t_0=20;

c_real=exp(-t_0/theta)*exp(-(t_0/eta)^m)*normcdf((mu-log(t_0))/sigma);%��ֵ
C_real=[exp(-t_0/theta);exp(-(t_0/eta)^m);normcdf((mu-log(t_0))/sigma)];
Result=zeros(1,1000);
Result_initial=zeros(1,1000);
Result_r=zeros(3,1000);
for i=1:1000
    i
    %����һ������
    data_1=-theta*log(rand(1,n_1));%ָ��
    data_2=eta*(-log(rand(1,n_2))).^(1/m);%������
    data_3=exp(normrnd(mu,sigma,1,n_3));%������̬
    data=[data_1;data_2;data_3];
    data_2=log(data_2);
    data_3=log(data_3);
    r_hat_1=exp(-t_0/(mean(data_1)));
    r_hat_2=exp(-exp((log(t_0)-mean(data_2))/(6^0.5*((n_2-1)*var(data_2)/n_2)^0.5/pi)-psi(1)));
    r_hat_3=normcdf((mean(data_3)-log(t_0))/(((n_3-1)*var(data_3)/n_3)^0.5));
    c=r_hat_1*r_hat_2*r_hat_3;%ϵͳ�ɿ��ԵĹ�����

    %��һ����ʼֵ
    r=[0.5,0.5,0.5];
    r_u=[1,1,1];
    r_l=[0,0,0];
    while  r_u-r_l>0.0001
        if G(r(1),r(2),r(3),n_1,n_2,n_3,c,aa)<0
            r_l=r;
            r=0.5*(r+r_u);
        end
        if G(r(1),r(2),r(3),n_1,n_2,n_3,c,aa)>0
            r_u=r;
            r=0.5*(r+r_l);
        end
        if G(r(1),r(2),r(3),n_1,n_2,n_3,c,aa)==0
            break
        end

    end
    Result_initial(i)=r(1)^3;
    %��һ����ʼֵ

    %�㷨��������
    beta=0.05;
    lambda=1;
    alpha=0.4;%��ʼ����
    %�㷨��������

    %Lagrange Newton�㷨
    P=1;
    P_temp=0;
    %while abs(P-P_temp)>0.00001
     k=1;
    for k=1:30%ֻ����һ������
    P_temp=P;
    P=(r(2)*r(3)-lambda*f_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N))^2+(r(1)*r(3)-lambda*f_2(r(1),r(2),r(3),n_1,n_2,n_3,c,N))^2+(r(1)*r(2)-lambda*f_3(r(1),r(2),r(3),n_1,n_2,n_3,c,N))^2+G(r(1),r(2),r(3),n_1,n_2,n_3,c,aa)^2;
    f12=(f_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N)-f_1(r(1),r(2)-0.01,r(3),n_1,n_2,n_3,c,N))/0.01;%����ƫ�������ò��̴���
    f13=(f_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N)-f_1(r(1),r(2),r(3)-0.01,n_1,n_2,n_3,c,N))/0.01;
    f23=(f_2(r(1),r(2),r(3),n_1,n_2,n_3,c,N)-f_2(r(1),r(2),r(3)-0.01,n_1,n_2,n_3,c,N))/0.01;
    W=[-lambda*df_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N),r(3)-lambda*f12,r(2)-lambda*f13,-f_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N);r(3)-lambda*f12,-lambda*df_2(r(1),r(2),r(3),n_1,n_2,n_3,c,N),r(1)-lambda*f23,-f_2(r(1),r(2),r(3),n_1,n_2,n_3,c,N);r(2)-lambda*f13,r(1)-lambda*f23,-lambda*df_3(r(1),r(2),r(3),n_1,n_2,n_3,c,N),-f_3(r(1),r(2),r(3),n_1,n_2,n_3,c,N);-f_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N),-f_2(r(1),r(2),r(3),n_1,n_2,n_3,c,N),-f_3(r(1),r(2),r(3),n_1,n_2,n_3,c,N),0];           
    y=-[r(2)*r(3)-lambda*f_1(r(1),r(2),r(3),n_1,n_2,n_3,c,N);r(1)*r(3)-lambda*f_2(r(1),r(2),r(3),n_1,n_2,n_3,c,N);r(1)*r(2)-lambda*f_3(r(1),r(2),r(3),n_1,n_2,n_3,c,N);-G(r(1),r(2),r(3),n_1,n_2,n_3,c,aa)];
    x=W\y;
    lambda_1=lambda-alpha*x(4);%%
    r_temp=r-alpha*x(1:3)';%%
    P_1=(r_temp(2)*r_temp(3)-lambda_1*f_1(r_temp(1),r_temp(2),r_temp(3),n_1,n_2,n_3,c,N))^2+(r_temp(1)*r_temp(3)-lambda_1*f_2(r_temp(1),r_temp(2),r_temp(3),n_1,n_2,n_3,c,N))^2+(r_temp(1)*r_temp(2)-lambda_1*f_3(r_temp(1),r_temp(2),r_temp(3),n_1,n_2,n_3,c,N))^2+G(r_temp(1),r_temp(2),r_temp(3),n_1,n_2,n_3,c,aa)^2;
    
    if P_1<=(1-beta*alpha)*P
        r=r_temp;
        lambda=lambda_1;
    else
        alpha=alpha/2;
    end
    
    end
%Lagrange Newton�㷨

     r_u=[1,1,1];
     r_l=[0,0,0];
    while  r_u-r_l>0.0001
       if G(r(1),r(2),r(3),n_1,n_2,n_3,c,aa)<0
            r_l=r;
            r=0.5*(r+r_u);
        end
        if G(r(1),r(2),r(3),n_1,n_2,n_3,c,aa)>0
            r_u=r;
            r=0.5*(r+r_l);
        end
        if G(r(1),r(2),r(3),n_1,n_2,n_3,c,aa)==0
            break
        end
    end
    R_L=r(1)*r(2)*r(3);
    GG(i)=G(r(1),r(2),r(3),n_1,n_2,n_3,c,aa);
    %ͳ�ƽ��
    Result(i)=R_L;
    Result_r(:,i)=r';
end
