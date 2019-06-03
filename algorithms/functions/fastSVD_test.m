%% test fast SVD

%generate data sets
n = 500;
p = 100;
R = 20;
X = randn(n,p);

% standard SVD
tic
for i = 1:10,
    [U,D,V] = svds(X,R);
end
svd_time = toc;

% fast SVD
tic
for i = 1:10,
    [U2,D2,V2] = fastSVD(X,R);
end
fastsvd_time = toc;

% test the results
norm(abs(U)-abs(U2),'fro')^2/norm(abs(U),'fro')^2
norm(abs(V)-abs(V2),'fro')^2/norm(abs(V),'fro')^2
norm(abs(D)-abs(D2),'fro')^2/norm(abs(D),'fro')^2

% compare the time 
svd_time
fastsvd_time


