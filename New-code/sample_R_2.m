function y=sample_R_2(r,n,N)%生成威布尔模型可靠度估计量的分布，这里的r是真值
global randnum_2
temp=randnum_2;
temp=log(log(temp.^(-1)));%标准极值分布样本
temp_mean=mean(temp);%W也即Q_1
temp_var=((n-1)/n)*var(temp);%V^2也即Q_2
y=exp(-exp((log(log(1/r))-temp_mean).*(1./temp_var)*pi/(6^0.5)-psi(1)));
end