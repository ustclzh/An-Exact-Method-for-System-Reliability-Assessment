function y=sample_R_1(r,n,N)%����ָ��ģ�Ϳɿ��ȹ������ķֲ��������r����ֵ
global randnum_1
temp=randnum_1;
temp=sum(temp.^2);%���ɶ�Ϊ2n_1�Ŀ����ֲ�
y=r.^(2*n./temp);
end
