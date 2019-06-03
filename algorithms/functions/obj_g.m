function [ out ] = obj_g( X2, A, B2 )
%least square loss for continuous data
out = 0.5*(norm(X2-A*B2','fro')^2);

end

