function [ out ] = obj_g_offset( X2, muX2, A, B2 )
%least square loss for continuous data
[m ~] = size(X2);
tmp = X2-A*B2'-ones(m,1)*muX2;
out = 0.5*(norm(tmp(~isnan(tmp)),'fro')^2);

end

