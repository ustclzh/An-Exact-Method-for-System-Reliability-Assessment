function y=sample_R_3(r,n,N)%����ָ��ģ�Ϳɿ��ȹ������ķֲ��������r����ֵ
global randnum_3
temp=randnum_3;
temp=sum(temp.^2);%���ɶ�Ϊ2n_1�Ŀ����ֲ�
y=r.^(2*n./temp);
end
