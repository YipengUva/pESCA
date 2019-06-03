function [cvErrors_mat,inits,outs] = ...
                        ESCA_modelSelection_KCV(dataSets,dataTypes,K,lambdas_md,fun,opts)
% Cross validation (CV) error based model selection procedure for the ESCA model with a group concave 
% penalty on the loading matrix to induce structured sparsity on the loading matrix to seperate 
% the global, local common and distinct variation in multiple mixed types data sets.
%
% Input:
%      dataSets: a cell array to hold the L data sets;
%      dataTypes: a cell array to hold the data types of the L data sets;
%      K: K-folds CV; 
%      lambdas_md: a vector of values of the tuning parameter lambda;
%      fun: the name of the used concave penalty;
%      opts.
%           tol_obj: tolerance for relative change of hist_obj, default:1e-6;
%           maxit: max number of iterations, default: 1000;
%           gamma:  hyper-parameter of penalty g(), default: 1;
%           R: the number of PCs, default: 50;
%           randomstart: initilization method, random (1), SCA(0), default: 0;
%           alphas: dispersion parameters of exponential dispersion families, default: 1. 
%           threPath: the option to generate thresholding path, default: 0;
%           type: 'L2', concave L2-norm penalty, only generate group sparsity;
%                 'L1', concave L1-norm penlaty, generate group sparsity, and slight elementwise sparsity;
%                 'biLevel', composite concave penalty, generate group and elementwise sparsity.
%
% Output:
%       cvErrors_mat: a length(lambdas_md)*(nDataSets+1)*K array to hold the CV errors;
%       inits: a length(lambdas_md)*1 cell structure to hold the possible initilizations for the final model;
%       outs: a length(lambdas_md)*1 cell structure to hold the output out during the model selection;	
					
% number of data sets, size of each data set
n_tries = length(lambdas_md); % number of values of lambda
nDataSets = length(dataSets); % number of data sets
n = nan(1,nDataSets);  % number of samples
d = nan(1,nDataSets);  % numbers of variables in different data sets
for i = 1:nDataSets, 
    [n(i), d(i)] = size(dataSets{i}); 
end
n    = unique(n);

% default dispersion parameters alphas
if isfield(opts, 'alphas'), alphas = opts.alphas; else alphas = ones(1,nDataSets); end

% create zero matrix to hold results
cvErrors_mat = zeros(n_tries, nDataSets+1, K); % +1 is used for the sum of all the Xi

% K folds CV
for k = 1:K
    opts_inner = opts;
    
    % split data sets into training set and test set
    ratio_mis = 0.1; % proportion to sample the test set  
    trainSets = cell(1,nDataSets); % training set
    testSets  = cell(1,nDataSets); % test set 
    indexSets = cell(1,nDataSets); % index of the test set
    for i = 1:nDataSets,
        
		% index out the i-th data set
        Xi = dataSets{i}; 
        dataType_Xi = dataTypes{i};
        
        % generate the index of the test set
		full_ind_vec = reshape(1:n*d(i),[n*d(i),1]);
        % if it is binary data, using hierachical sampling
        if strcmp(dataType_Xi, 'Bernoulli'),
            ones_ind_vec  = full_ind_vec(Xi == 1);
            zeros_ind_vec = full_ind_vec(Xi == 0);
            index_Xi_ones  = randsample(ones_ind_vec, ceil(ratio_mis*length(ones_ind_vec)));
            index_Xi_zeros = randsample(zeros_ind_vec, ceil(ratio_mis*length(zeros_ind_vec)));
  
            % test the sampled samples
            if or(not(all(Xi(index_Xi_ones)==1)),not(all(Xi(index_Xi_zeros)==0))),
                disp('the hierachical sampling does not work');
            end
            
            % combine the index_Xi
            index_Xi = [index_Xi_ones;index_Xi_zeros];
            
        else % if nonBinary data, using randomly sampling
            non_NaN_mat  = 1-isnan(Xi);
            non_NaN_vec  = non_NaN_mat(:);
            non_NaN_ind_vec = full_ind_vec(non_NaN_vec > 0);
            index_Xi = randsample(non_NaN_ind_vec, ceil(ratio_mis*length(non_NaN_ind_vec)));
        end
        
        % generate the train set
        Xi_train = Xi; Xi_train(index_Xi) = nan;
        trainSets{i} = Xi_train;
        
        % generate the test set
        Xi_test = Xi(index_Xi);
        
        testSets{i}  = Xi_test;
        indexSets{i} = index_Xi;
    end
    
    % save the parameters during the model selection
    if (k==K)
        inits = cell(n_tries,1);
		outs  = cell(n_tries,1);
    end
        
    % model selection
    for j = 1:n_tries
        lambda = lambdas_md(j);
	
        % using the training set to construct a ESCA model
        lambdas = lambda*ones(1,nDataSets);
	    [mu,A,B,~,out_inner] = ESCA_group_concave(trainSets,dataTypes,lambdas,fun,opts_inner);
	    if (out_inner.iter <= 2),
            disp('less than 3 iteration is used.');
        end
        
	    % warm start
        opts_inner.mu0 = mu;
        opts_inner.A0  = A;
        opts_inner.B0  = B;
        
        if (k==K),
            inits{j,1} = opts_inner;
			outs{j,1}  = out_inner;
        end

	    % compute the test error
        ThetaHat  = ones(n,1)*mu' + A*B';
        testError_vec = zeros(1,nDataSets);
        
        for i = 1:nDataSets,
            % index out ThetaHat_Xi
            columns_Xi = (sum(d(1:(i-1)))+1):sum(d(1:i));
            ThetaHat_Xi = ThetaHat(:,columns_Xi);
            
            % compute the CV error
            index_Xi = indexSets{i};
            Xi_test  = testSets{i};
            dataType_Xi = dataTypes{i};
            if strcmp(dataType_Xi, 'Gaussian'),
                testError_Xi = (1/alphas(i))*0.5*norm(Xi_test - ThetaHat_Xi(index_Xi),'fro')^2;
            elseif strcmp(dataType_Xi, 'Bernoulli'),
                testError_Xi = (1/alphas(i))*obj_f_logistic(Xi_test, ThetaHat_Xi(index_Xi));
            end
            testError_vec(1,i) = testError_Xi;
        end
        
        cvErrors_tmp  = [sum(testError_vec),testError_vec];
        cvErrors_mat(j,:,k) = cvErrors_tmp; 
    end
end

end

