function [cvErrors] = SVD_CV_modelSelection(X,K,Rs,opts,missingRatio)
% Missing value based CV model selection approach.
% Randomly select 10% percent elements as missing values.
% After that a EM-SVD model is constructed to estimate the prediction error.
%
% Input:
%      X,opts: same as SVD_missingValues.m
%      K: K times estimation
%      Rs: a vector R (the number of components)
%
% Output:
%       cvErrors: length(Rs)*K vector to hold cv errors 

% default parameter
if(nargin<5), missingRatio = 0.1; end

% structure to hold results
length_Rs = length(Rs);
cvErrors  = nan(length_Rs,K);

% number of non-missing elements
[m,n] = size(X);
mn = m*n;
mn_nonNaN = sum(sum(1-isnan(X)));

for k = 1:K,
    % seperate the Xtest and Xtrain
    % taken into account the potential problem of NaN
    full_ind_vec = reshape(1:mn,[mn,1]);
    non_NaN_mat  = 1-isnan(X);
    non_NaN_vec  = non_NaN_mat(:);
    non_NaN_ind_vec = full_ind_vec(non_NaN_vec > 0);
    X_index = randsample(non_NaN_ind_vec, ceil(missingRatio*mn_nonNaN));
    X_train = X; 
    X_train(X_index) = nan;
    X_test  = X(X_index);
   
    % use opts_inner to do warm start
    opts_inner = opts;

    % for loop
    for j = length_Rs:-1:1
        R = Rs(j);
	
	    % using the remaining data to construct a SVD model
        [U,S,V,~] = SVD_missingValues(X_train,R,opts_inner);
	    ZHat = U*S*V'; 
        
        % warm start
        opts_inner.U0 = U;
        opts_inner.S0 = S;
		opts_inner.V0 = V;
            
	    % extract the estimated parameters for the prediction of missing elements
        X_pred = ZHat(X_index);
	
	    % compute the prediction error
	    cvErrors(j,k) = 0.5*norm(X_test-X_pred,'fro')^2;
    end
end

end

