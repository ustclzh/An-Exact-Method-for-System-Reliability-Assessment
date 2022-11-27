



clear;
clc;
n_1=10;
n_2=10;
n_3=10;
N=100000;


%所有随机数在这里生成
global randnum_1 randnum_2 randnum_31 randnum_32 Q_1 Q_2 Q_31 Q_32
randnum_1=normrnd(0,1,n_1,N);
randnum_2=rand(n_2,N);
randnum_31=normrnd(0,1,1,N);
randnum_32=normrnd(0,1,n_3-1,N);
Q_1=normrnd(0,1,n_1,N);
Q_2=rand(n_2,N);
Q_31=normrnd(0,1,1,N);
Q_32=normrnd(0,1,n_3-1,N);


aa=0.2; %置信水平alpha
eta=200;
m=2;
theta=200;
sigma=2;
mu=3;
t_0=20;
c_real=exp(-t_0/theta)*exp(-(t_0/eta)^m)*normcdf((mu-log(t_0))/sigma);%真值

for i=1:1000

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


beta=0.05;
r=[0.5,0.5,0.5];
r_u=[1,1,1];
r_l=[0,0,0];
while        r_u-r_l>0.000001
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

RR_L(i)=r(1)^3;
end




