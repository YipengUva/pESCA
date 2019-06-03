function [varExp_total,varExp_PCs] = varExp_Gaussian(X,mu,A,B,Q)
% compute the variation expalined of Gaussian distributed data
% input:
%      X: quantitative data set
%      mu: column offset term
%      A: score matrix
%      B: loading matrix
%
% output:
%      varExp_total: variation explained.
%      varExp_PCs: variation explained for each PCs.

% parameters used 
[m,~]       = size(X);

% compute the loglikelihood of mle and null model
X_centered  = X - ones(m,1)*mu';
QX          = Q.*X_centered;
QXt         = QX';

% likelihood of null model
l_null  = norm(QX,'fro')^2; % null model 

% likelihood of full model
E_hat   = X_centered - A*B';
QE_hat  = Q.*E_hat;
l_model = norm(QE_hat, 'fro')^2; % full model

% compute the least squares of an individual PC
[~,R]  = size(B);
l_PCs  = zeros(1,R);

for r=1:R
    Ar = A(:,r); Br = B(:,r); Brt = Br';
    QPCr = Q.*(Ar*Brt);
	l_PCs(1,r) = l_null - 2*((Brt*QXt)*Ar) + Ar'*(QPCr*Br);
end

% compute variation explained by each PC
varExp_PCs = (1 - l_PCs./l_null)*100;

% total variation explained
varExp_total = (1 - l_model/l_null)*100;

end

