function y=G(r_1,r_2,r_3,n_1,n_2,n_3,c,alpha)%生成R的估计量的分布，用于计算G的近似分布。这里的r_1,r_2,r_3都是真值
%N=100000;
global N;
temp=sample_R_1(r_1,n_1,N).*sample_R_2(r_2,n_2,N).*sample_R_3(r_3,n_3,N);%串联系统
y=sum(temp>c)/N-alpha;
end


