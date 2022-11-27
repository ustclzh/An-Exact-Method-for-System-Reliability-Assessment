function y=sample_R_3(r,n,N)%生成对数正态模型可靠度估计量的分布，这里的r是真值
global randnum_31 randnum_32 
temp_mean=randnum_31;%标准正态
temp=randnum_32 ;
temp_var=sum(temp.^2);%自由度为n-1的卡方分布
r_1=norminv(r,0,1);%求分位点
y=normcdf(((r_1*n^0.5+temp_mean)./(temp_var.^0.5)),0,1);
end