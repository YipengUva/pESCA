function [U,D,V] = fastSVD(X,R)
% fast SVD use the kernel trick for the SVD of a fat matrix
% used for the matrix X, p >> n

% size of the matrix 
[n,p] = size(X);
Xt = X';

% select the svd method
if (n < p),
    XXt = X*Xt;
    [U,Dsquare] = eigs(XXt,R);
    D_vec = sqrt(diag(Dsquare));
    D = diag(D_vec);
    V = (Xt*U)*diag(1./D_vec);
else
    [U,D,V] = svds(X,R);
end

end

