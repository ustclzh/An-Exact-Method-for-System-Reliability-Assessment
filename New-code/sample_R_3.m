function y=sample_R_3(r,n,N)%���ɶ�����̬ģ�Ϳɿ��ȹ������ķֲ��������r����ֵ
global randnum_31 randnum_32 
temp_mean=randnum_31;%��׼��̬
temp=randnum_32 ;
temp_var=sum(temp.^2);%���ɶ�Ϊn-1�Ŀ����ֲ�
r_1=norminv(r,0,1);%���λ��
y=normcdf(((r_1*n^0.5+temp_mean)./(temp_var.^0.5)),0,1);
end