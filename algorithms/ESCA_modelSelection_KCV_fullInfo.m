function [cvErrors_mat,ranks_mat,RVs_mat,RMSEs_mat,inits,outs] = ...
                        ESCA_modelSelection_KCV_fullInfo(ds,dataTypes,K,lambdas_md,dataSimulation,fun,opts)

% CV error based model selection procedure when the simulated parameters 
% are avaliable.
%
% Input:
%      ds: the number of variables of each data set;
%      dataSimulation: the output structure of the data simulation process;
%      dataTypes,K,lambdas_md,fun,opts: same as the ESCA_modelSelection_KCV.m;
%
% Output:
%       ranks_mat: a length(lambdas_md)*7 matrix to hold the rank estimations 
%                  for the global common (C123), local common (C12,C13,C23) 
%                  and distinct (D1,D2,D3) structures;
%       RVs_mat: a length(lambdas_md)*7 matrix to hold the RV coefficients 
%                for the global common (C123), local common (C12,C13,C23) 
%                and distinct (D1,D2,D3) structures;
%       RMSEs_mat: a length(lambdas_md)*5 matrix to hold the RMSEs in 
%                  estimating \Theta, \{\Theta_l \}_1^3, \mu; 
%       cvErrors_mat,inits,outs: same as the ESCA_modelSelection_KCV.m.      

% default parameters
if isfield(opts, 'alphas'), alphas = opts.alphas; else alphas = ones(1,length(dataTypes)); end
n_tries = length(lambdas_md);

% index out the simulated data sets 
X      = dataSimulation.X;
[n,~]  = size(X);
X1 = X(:,1:ds(1));
X2 = X(:,(ds(1)+1):(sum(ds(1:2))));
X3 = X(:,(sum(ds(1:2))+1):end);
dataSets  = {X1, X2, X3};
nDataSets = length(dataSets);

% create zero matrix to hold results
cvErrors_mat = zeros(n_tries, 4, K);
ranks_mat    = zeros(n_tries, 7, K);
RVs_mat      = zeros(n_tries, 7, K);
RMSEs_mat    = zeros(n_tries, 5, K);

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
		full_ind_vec = reshape(1:n*ds(i),[n*ds(i),1]);
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

    if (k==K)
        inits = cell(n_tries,1);
		outs  = cell(n_tries,1);		
    end
        
    % model selection
    for j = 1:n_tries
        lambda = lambdas_md(j);
	
        % using the training set to construct a ESCA model
        lambdas = lambda*ones(1,nDataSets);
	    [mu,A,B,S,out_inner] = ESCA_group_concave(trainSets,dataTypes,lambdas,fun,opts_inner);
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
        
		% model evulation with respect to the simulated parameters
		[RVs_structures,Ranks_structures,RMSEs_parameters] = ...
		      ESCA_evaluation_metrics(mu,A,B,S,ds,dataSimulation);

        RVs_mat(j,:,k)   = RVs_structures;
    	ranks_mat(j,:,k) = Ranks_structures;
	    RMSEs_mat(j,:,k) = RMSEs_parameters;
		
	    % compute the test error
        ThetaHat  = ones(n,1)*mu' + A*B';
        testError_vec = zeros(1,nDataSets);
        
        for i = 1:nDataSets,
            % index out ThetaHat_Xi
            columns_Xi = (sum(ds(1:(i-1)))+1):sum(ds(1:i));
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

