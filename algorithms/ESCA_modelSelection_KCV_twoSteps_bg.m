function [cvErrors_g,cvErrors_b,lambdas_opt,opts_init] = ...
                        ESCA_modelSelection_KCV_twoSteps_bg(dataSets,dataTypes,...
                        lambda_g0,lambdas_md_g,lambdas_md_b,fun,opts)
% CV error based model selection procedure for tuning the model on Gaussian and Bernoulli mixture data sets. 
% First select lambda_b for Bernoulli data sets then select lambda_g for Gaussian data sets.
%
% Input:
%      lambda_g0: initial value of lambda_g for Gaussian data sets;
%      lambdas_md_g: vaules of lambda_g for Gaussian data sets;
%      lambdas_md_b: vaules of lambda_b for Bernoulli data sets;
%      dataSets,dataTypes,fun,opts: same as ESCA_modelSelection_KCV.m. 
%
% Output:
%       cvErrors_g: a length(lambdas_md)*nGuassian matrix to hold the CV errors for tuning lambda_g;
%       cvErrors_b: a length(lambdas_md)*nBernoulli matrix to hold the CV errors for tuning lambda_b;
%       lambdas_opt: optimal values of lambdas for a ESCA model;
%       opts_init: the initilization of the selected model.

% number of data sets, size of each data set
nDataSets = length(dataSets); % number of data sets
n = nan(1,nDataSets);  % number of samples
d = nan(1,nDataSets);  % numbers of variables in different data sets
for i = 1:nDataSets, 
    [n(i), d(i)] = size(dataSets{i}); 
end
n    = unique(n);

% index out the how many Gaussian data sets and how many number nonGuassian data sets
length_g = 0;
for i = 1:nDataSets,
    dataType_Xi = dataTypes{i};
	if strcmp(dataType_Xi, 'Gaussian'),
	    length_g = length_g + 1;
	end
end
length_b = nDataSets - length_g;

% create zero matrix to hold results
lambdas_opt    = zeros(1,nDataSets);

% model selection using CV
opts_inner = opts;
    
% split data sets into trainSets and testSets
ratio_mis = 0.1;
trainSets = cell(1,nDataSets);
testSets  = cell(1,nDataSets);
indexSets = cell(1,nDataSets);
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
    
% model selection for nonGuassian data sets
% first fix lambda_g0, optimize lambda_b
lambda_g0 = ones(1,length_g)*lambda_g0;
[cvErrors_b,inits_b] = ...
                    ESCA_modelSelection_KCV_twoSteps_nonGaussian(trainSets,testSets,indexSets,dataTypes,...
                    lambda_g0,lambdas_md_b,fun,opts_inner);	
						
% select the optimal value of lambda_b for nonGuassian data sets
nonGaussian_cvErrors     = cvErrors_b(:,(length_g+1):end);
sum_nonGaussian_cvErrors = sum(nonGaussian_cvErrors,2);
[~,index_b]  = min(sum_nonGaussian_cvErrors);
lambda_b_opt = ones(1,length_b)*lambdas_md_b(index_b);
opts_inner   = inits_b{index_b,1};
    
% model selection for Gaussian data sets
lambda_b0 = lambda_b_opt;
[cvErrors_g,inits_g] = ...
                    ESCA_modelSelection_KCV_twoSteps_Gaussian(trainSets,testSets,indexSets,dataTypes,...
                    lambda_b0,lambdas_md_g,fun,opts_inner);
	
% select the optimal value of lambda_b for nonGuassian data sets
Gaussian_cvErrors     = cvErrors_g(:,1:length_g);
sum_Gaussian_cvErrors = sum(Gaussian_cvErrors,2);
[~,index_g]  = min(sum_Gaussian_cvErrors);
lambda_g_opt = ones(1,length_g)*lambdas_md_g(index_g);
opts_init    = inits_g{index_g,1};
	
% save the results
lambdas_opt(1,:)  = [lambda_g_opt, lambda_b_opt];

end



