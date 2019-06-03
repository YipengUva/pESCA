function [cvErrors_g,inits_g] = ...
                        ESCA_modelSelection_KCV_twoSteps_Gaussian(trainSets,testSets,indexSets,...
						dataTypes,lambda_b0,lambdas_md_g,fun,opts_inner)
% used in ESCA_modelSelection_KCV_twoSteps_bg.m to turn lambda_g for Gaussian data sets.	
			
% number of data sets, size of each data set
nDataSets = length(trainSets); % number of data sets
n = nan(1,nDataSets);  % number of samples
d = nan(1,nDataSets);  % numbers of variables in different data sets
for i = 1:nDataSets, 
    [n(i), d(i)] = size(trainSets{i}); 
end
n    = unique(n);
if isfield(opts_inner, 'alphas'), alphas = opts_inner.alphas; else alphas = ones(1,nDataSets); end

% index out the how many Gaussian data sets and how many number nonGuassian data sets
length_g = 0;
for i = 1:nDataSets,
    dataType_Xi = dataTypes{i};
	if strcmp(dataType_Xi, 'Gaussian'),
	    length_g = length_g + 1;
	end
end

% model selection of nonGaussian data sets
% create zero matrix to hold results
n_tries_g  = length(lambdas_md_g);
inits_g    = cell(n_tries_g,1);
cvErrors_g = zeros(n_tries_g, nDataSets);

for j = 1:n_tries_g,
    lambda_g = ones(1,length_g)*lambdas_md_g(j);
	
	lambdas = [lambda_g, lambda_b0];
	[mu,A,B,~,out_inner] = ESCA_group_concave(trainSets,dataTypes,lambdas,fun,opts_inner);
	if (out_inner.iter <= 2),
        disp('less than 3 iteration is used.');
    end
						
	% warm start
    opts_inner.mu0 = mu;
    opts_inner.A0  = A;
    opts_inner.B0  = B;
	
	inits_g{j,1} = opts_inner;
	
	% compute the test error
    ThetaHat  = ones(n,1)*mu' + A*B';
    testError_vec = zeros(1,nDataSets);
        
    for i = 1:nDataSets,
        % index out ThetaHat_Xi
        columns_Xi  = (sum(d(1:(i-1)))+1):sum(d(1:i));
        ThetaHat_Xi = ThetaHat(:,columns_Xi);
        alpha_Xi    = alphas(i);
        
        % compute the CV error
        index_Xi = indexSets{i};
        Xi_test  = testSets{i};
        dataType_Xi = dataTypes{i};
        if strcmp(dataType_Xi, 'Gaussian'),
            testError_Xi = (1/alpha_Xi)*0.5*norm(Xi_test - ThetaHat_Xi(index_Xi),'fro')^2;
        elseif strcmp(dataType_Xi, 'Bernoulli'),
            testError_Xi = (1/alpha_Xi)*obj_f_logistic(Xi_test, ThetaHat_Xi(index_Xi));
        end
        testError_vec(1,i) = testError_Xi;
    end

    cvErrors_g(j,:) = testError_vec; 
end

end

