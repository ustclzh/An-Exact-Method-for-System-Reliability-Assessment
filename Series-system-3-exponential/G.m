function y=G(r_1,r_2,r_3,n_1,n_2,n_3,c,alpha)%����R�Ĺ������ķֲ������ڼ���G�Ľ��Ʒֲ��������r_1,r_2,r_3������ֵ
%N=100000;
global N;
temp=sample_R_1(r_1,n_1,N).*sample_R_2(r_2,n_2,N).*sample_R_3(r_3,n_3,N);%����ϵͳ
y=sum(temp>c)/N-alpha;
end


