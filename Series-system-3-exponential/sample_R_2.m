function y=sample_R_2(r,n,N)%����ָ��ģ�Ϳɿ��ȹ������ķֲ��������r����ֵ
global randnum_2
temp=randnum_2;
temp=sum(temp.^2);%���ɶ�Ϊ2n_1�Ŀ����ֲ�
y=r.^(2*n./temp);
end
