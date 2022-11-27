function y=R_weibull(r,Q_1,Q_2)%模拟威布尔分布R的分布，参数总共2个，样本量n以及真值r，N是模拟的总数。n和N分别在矩阵Q_1,Q_2中体现

y=exp(-exp((log(log(1/r))-Q_1)./(Q_2)*pi/(6^0.5)-psi(1)));
end