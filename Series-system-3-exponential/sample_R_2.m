function y=sample_R_2(r,n,N)%生成指数模型可靠度估计量的分布，这里的r是真值
global randnum_2
temp=randnum_2;
temp=sum(temp.^2);%自由度为2n_1的卡方分布
y=r.^(2*n./temp);
end
