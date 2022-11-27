


function [c,ceq]=confun(x)
c=[];
global c_1 aa;
ceq=G(x(1),x(2),x(3),10,10,10,c_1,aa);
end