function y=sample_R_2(r,n,N)%����������ģ�Ϳɿ��ȹ������ķֲ��������r����ֵ
global randnum_2
temp=randnum_2;
temp=log(log(temp.^(-1)));%��׼��ֵ�ֲ�����
temp_mean=mean(temp);%WҲ��Q_1
temp_var=((n-1)/n)*var(temp);%V^2Ҳ��Q_2
y=exp(-exp((log(log(1/r))-temp_mean).*(1./temp_var)*pi/(6^0.5)-psi(1)));
end