function y=f_2(r_1,r_2,r_3,n_1,n_2,n_3,c,N)%输出关于r_1的密度估计
 %R_1=sample_R_1(r_1,n_1,N);
 R_1=sample_R_1(r_1,n_1,N);
 R_3=sample_R_3(r_3,n_3,N);
 
temp=c./(R_1.*R_3);
global Q_2;
temp1=Q_2;
Q=sum(temp1.^2); %卡方分布

%所用核函数为（3/2-1/2*u^2)*exp(-u^2/2)/((2*pi)^0.5)

X=temp.^(Q/(2*n_2));%模拟数据，用于密度估计
h=0.01;
y=1/(N*h)*sum((3/2-1/2*(X-r_2).^2).*exp(-(X-r_2).^2/2)/((2*pi)^0.5));
end