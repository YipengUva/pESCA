function y = GDP_sg(x,gamma,lambda)
% supergradient of GDP penalty
% x>=0

x = abs(x) ;
y = lambda./(gamma + x);
