function [result] = fastTrace(X,Y)
% fast trace function
% if n>p, trace(XY) = trace(YX);
% if n<p, trace(XY) = trace(XY);

% size of the matrix 
[n,p] = size(X); %Y p*n

if (n > p),
    result = trace(Y*X);
else
    result = trace(X*Y);
end

end

