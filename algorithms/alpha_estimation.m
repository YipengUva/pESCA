function [est_mean,est_std,R_CV] = alpha_estimation(X,K,Rs,opts)
% alpha estimation procedure using the PCA model.
% The rank of the PCA model is selected based on a CV error based model 
% selection procedure.
%
% Input:
%      X: a quantitative data set;
%      K: K-fold CV;
%      Rs: the searching range of the number of components;
%      opts.
%           tol_obj: tolerance for relative change of hist_obj, default:1e-4;
%           maxit: max number of iterations, default: 1000;
%
% Output:
%       est_mean: mean of the K times alpha estimations; 
%       est_std: standard deviation of K times alpha estimations;
%       R_CV: selected ranks in K times alpha estimations;

if(nargin<4), opts = []; end

% check if the whole row or whole column is missing
% index out the rows and columns not fully missing
W = 1-isnan(X);
X = X(sum(W,2)>0,sum(W,1)>0);

% first center the data sets
[m,n]  = size(X);
X_mean = mean(X,1,'omitnan');
X      = X - ones(m,1)*X_mean;

% parameters
W = 1-isnan(X);
mn_nonNaN = sum(sum(W));

% model selection
[cvErrors] = SVD_CV_modelSelection(X,K,Rs,opts);

% select the number of components
[~,CV_index] = min(cvErrors,[],1);
R_CV = Rs(CV_index);

% estimate the noise level
alphas = nan(1,length(R_CV));

% first do a SVD on X to accelerate computation
if all(all(1-isnan(X))),
   [U,S,V] = fastSVD(X,max(R_CV));
else
    opts.tol_obj = 1e-6;
    [U,S,V,~] = SVD_missingValues(X,max(R_CV),opts);
end

for i=1:length(R_CV),
    R_CV_tmp = R_CV(i);
    Z_hat    = U(:,1:R_CV_tmp)*S(1:R_CV_tmp,1:R_CV_tmp)*V(:,1:R_CV_tmp)';
    DF = mn_nonNaN - (m + n)*R_CV_tmp;
    X_nonNaN = X;
    X_nonNaN(isnan(X)) = 0;
    sigmaSqure = (1/DF)*norm(X_nonNaN-Z_hat,'fro')^2;
    alphas(1,i) = sigmaSqure;
end

est_mean = sum(alphas)/K;
est_std  = std(alphas);

% CV error plot
figure
plot(Rs, cvErrors, '-o');
title('CV errors');
xlabel('ranks')

end

