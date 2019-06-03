function y = GDP(x,gamma,lambda)
% smooth scad penalty for singular values
% x>=0

x = abs(x) ;
y = lambda*log(1+x./gamma);

