function [varExp_total, varExp_PCs] = GSCA_varExp(X1,X2,R,mu,Z,sigmaSquare)
% compute the variation expalined
% input:
%      X1: binary data set
%      X2: quantitative data set
%      R: number of PCs;
%      mu: column offset term
%      Z: low rank approximation
% output:
%      varExp_total: variation explained.
%      varExp_PCs: variation explained for each PCs.

% parameters used 
[m,n1]       = size(X1);
[~,n2]       = size(X2);
Q1  = 1-isnan(X1); % weighting matrix for accounting missing values

offsetX1 = ones(m,1)*mu(1:n1,:)';
offsetX2 = ones(m,1)*mu(n1+1:end,:)';
ThetaHat = ones(m,1)*mu' + Z;
Theta1Hat = ThetaHat(:,1:n1);
Theta2Hat = ThetaHat(:,n1+1:end);

% compute the log likelihood of mle and null model
l_model = obj_logistic(X1,Theta1Hat,Q1) +...
    (1/(2*sigmaSquare))*norm((X2-Theta2Hat),'fro')^2 ...
        + 0.5*m*n2*log(2*pi*sigmaSquare);
l_null = obj_logistic(X1,offsetX1,Q1) +...
    (1/(2*sigmaSquare))*norm((X2-offsetX2),'fro')^2 ...
        + 0.5*m*n2*log(2*pi*sigmaSquare);
    
% deviance ratio
varExp_total = 1 - l_model/l_null;

% compute the log likelihood of a model with an individual PC
[U,S,V] = svds(Z,R);
l_PCs  = zeros(1,R);

for r=1:R
    Theta_PC  = ones(m,1)*mu' + U(:,r)*S(r,r)*V(:,r)';
    Theta1_PC = Theta_PC(:,1:n1);
    Theta2_PC = Theta_PC(:,n1+1:end);
    l_PCs(1,r) = obj_logistic(X1,Theta1_PC,Q1) +...
    (1/(2*sigmaSquare))*norm((X2-Theta2_PC),'fro')^2 ...
        + 0.5*m*n2*log(2*pi*sigmaSquare);
end

% compute variation explained by each PC
varExp_PCs = 1 - l_PCs./l_null;

end

