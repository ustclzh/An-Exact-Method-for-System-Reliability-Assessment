function y=df_3(r_1,r_2,r_3,n_1,n_2,n_3,c,N)%输出关于r_1的密度估计
 R_1=sample_R_1(r_1,n_1,N);
 R_2=sample_R_3(r_2,n_2,N);
 %R_3=sample_R_3(r_3,n_3,N);
 
temp=c./(R_1.*R_2);
%alpha=sum(temp>=1)/N;
global Q_31 Q_32;
Q1=Q_31;
temp1=Q_32;
Q2=sum(temp1.^2);
%所用核函数为（3/2-1/2*u^2)*exp(-u^2/2)/((2*pi)^0.5)

X=normcdf(((norminv(temp,0,1).*(Q2).^0.5-Q1)/(n_3^0.5)),0,1);%模拟数据，用于密度估计
X(temp>=1)=2;

h=0.01;
y=1/(N*h)*sum((X-r_3).*exp(-(X-r_3).^2/2)/((2*pi)^0.5)-(3/2-1/2*(X-r_3).^2).*((X-r_3)/((2*pi)^0.5)).*exp(-(X-r_3).^2/2)/((2*pi)^0.5));
end