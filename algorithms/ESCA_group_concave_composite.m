function [mu,A,B,S,out] = ESCA_group_concave_composite(dataSets,dataTypes,lambdas,fun,opts)

% Suppose B_lr indicates the r-th column of the loading matrix B_l of the l-th data set
% X_l, and g() is the concave penalty. The composite concave penalty g(\sum_j g(|B_lrj|)) is
% used to induce group sparsity. 
% Furthermore, elements in B_lr are regularized in the same way as the concave penalty g().

% default parameters
if isfield(opts, 'tol_obj'),tol_obj   = opts.tol_obj;  else tol_obj  = 1e-6;  end
if isfield(opts, 'maxit'),   maxit    = opts.maxit;    else maxit    = 1000;  end
if isfield(opts, 'gamma'),   gamma    = opts.gamma;    else gamma    = 1;     end
if isfield(opts, 'R'),       R        = opts.R;        else R        = 50;    end 
if isfield(opts, 'random_start'), random_start = opts.random_start;...
else random_start=0; end
if isfield(opts, 'threPath'),threPath = opts.threPath;  else threPath= 0;  end
if isfield(opts, 'quiet'),    quiet   = opts.quiet;     else quiet   = 0;  end

% concave penalty function and its supergradient
hfun    = str2func(fun);         % penalty function 
hfun_sg = str2func([fun '_sg']); % super gradient of penalty function

% number of data sets, size of each data set
nDataSets = length(dataSets); % number of data sets
n = nan(1,nDataSets);  % number of samples
d = nan(1,nDataSets);  % numbers of variables in different data sets
for i = 1:nDataSets, 
    [n(i), d(i)] = size(dataSets{i}); 
end
sumd = sum(d); % total number of variables
n    = unique(n);

% form full data set X, X = [X1,...Xl,...XL]
% form full weighting matrix W, W = [W1,...Wl,...WL]
X = nan(n, sumd);
W = nan(n, sumd);
for i=1:nDataSets,
    columns_Xi = (sum(d(1:(i-1)))+1):sum(d(1:i));
	X_i = dataSets{i};
	W(:,columns_Xi)  = 1-isnan(X_i);
	X_i(isnan(X_i))  = 0;   % set missing values to be 0; then X_i =  W_i.*X_i;
	X(:,columns_Xi)  = X_i;
end

% default dispersion parameters alphas
if isfield(opts, 'alphas'), alphas = opts.alphas; else alphas = ones(1,nDataSets); end

% initialization
if(isfield(opts, 'A0'))     % use imputted initialization
    mu0 = opts.mu0; mu0 = mu0';
	A0  = opts.A0;
	B0  = opts.B0;
	R   = size(B0,2);
elseif (random_start == 1), % use random initialization
    mu0 = zeros(1,sumd);
    A0  = randn(n,R);
    B0  = randn(sumd,R);
elseif (random_start == 0), % use SCA model as initialization
    if (not(ismember('Poisson', dataTypes)))
        mu0 = mean(X,1);
        [U,D,V] = fastSVD(X-ones(n,1)*mu0, R);
        A0 = U;
        B0 = V*D;
    else
        X_tmp = X;
        for i=1:nDataSets, % log transformation applied to count data
          dataType = dataTypes{i};
          columns_Xi = (sum(d(1:(i-1)))+1):sum(d(1:i));
          if strcmp(dataType, 'Poisson'),
              X_tmp(:,columns_Xi) = log(X(:, columns_Xi)+1);
          end
        end
        mu0 = mean(X_tmp,1);
        [U,D,V] = fastSVD(X_tmp-ones(n,1)*mu0, R);
        A0 = U;
        B0 = V*D;
    end  
end
Theta0 = ones(n,1)*mu0 + A0*B0';

% initial value of loss function
f_obj0 = 0;
g_obj0 = 0;
Sigmas0 = nan(nDataSets,R); % sigma_{lr} = \sum g(blr)
for i = 1:nDataSets,
    dataType = dataTypes{i};
    columns_Xi = (sum(d(1:(i-1)))+1):sum(d(1:i));
	X_i        = X(:, columns_Xi);
	W_i        = W(:, columns_Xi);
	Theta0_i   = Theta0(:,columns_Xi);
	alpha_i    = alphas(i);
	weight_i   = d(i);  % weight when composited concave penalty is used
	lambda_i   = lambdas(i)*weight_i;
    
	log_partition = str2func([dataType '_b']); % specify log-partiton function
	W_it = W_i';
    Theta0_it = Theta0_i';
	f_obj0 = f_obj0 + (1/alpha_i)*...
	    (fastTrace(W_it,log_partition(Theta0_i)) - fastTrace(Theta0_it,X_i));
	
	for r = 1:R,
	    sigma0_ir    = sum(hfun(abs(B0(columns_Xi,r)),gamma,1)); % sigma_{lr} = \sum g(blrj)
	    Sigmas0(i,r) = sigma0_ir; 
		g_obj0 = g_obj0 + lambda_i*hfun(sigma0_ir,gamma,1);
    end	
end

obj0 = f_obj0 + g_obj0;   % objective + penalty
out.f_obj(1)    = f_obj0; % objective
out.g_obj(1)    = g_obj0; % penalty
out.hist_obj(1) = obj0;

% iterations
for k = 1:maxit
    if (quiet==0), fprintf('%d th iteration\n',k); end
    
	%--- form Hk ---
	%--- update mu ---
	%--- form JHk ---
	%--- form scaled_JHk ---
	%--- form scaled_Bk  ---
	JHk = nan(n,sumd);
	scaled_JHk = nan(n,sumd);
	scaled_Bk  = nan(sumd,R);
	mu   = nan(1,sumd);
	rhos = nan(1,nDataSets); % form rhos, the Lipshize constant for each data types
	for i = 1:nDataSets,
	    dataType   = dataTypes{i};
	    columns_Xi = (sum(d(1:(i-1)))+1):sum(d(1:i));
		X_i        = X(:,columns_Xi);
		Theta0_i   = Theta0(:, columns_Xi);
		W_i        = W(:,columns_Xi);
		log_partition_g = str2func([dataType '_bg']); % gradient of log partition function
		
		% form rhos, the Lipshize constant for each data types
		if strcmp(dataType, 'Bernoulli')
	        rho_i = 0.25;
	    elseif strcmp(dataType, 'Gaussian')
	        rho_i = 1;
	    elseif strcmp(dataType, 'Poisson')
	        rho_i = max(max(exp(Theta0_i)));
	    end
		rhos(i) = rho_i;
		
		% form Hk_i 
		Hk_i = Theta0_i - (1/rho_i)*(W_i .* (log_partition_g(Theta0_i)-X_i));
		
		% update mu_i
		mu_i = mean(Hk_i,1);
		mu(1,columns_Xi) = mu_i;
		
		% form JHk_i
		JHk_i = Hk_i - ones(n,1)*mu_i;
        JHk(:,columns_Xi) = JHk_i;
		
		% for scaled_JHk_i, scaled_Bk_i
		alpha_i = alphas(i);
		c_i     = sqrt(rho_i/alpha_i);
		
		scaled_JHk(:,columns_Xi) = c_i*JHk_i;
		scaled_Bk(columns_Xi,:)  = c_i*B0(columns_Xi,:);
	end
	
	%--- update A ---
	[U,~,V] = fastSVD(scaled_JHk*scaled_Bk, R);
	A = U*V';
	
	%--- update B ---
	B = nan(sumd,R);
	Sigmas    = nan(nDataSets,R);  % sigma_{lr} = \sum g(blrj)
    for i=1:nDataSets,
	    columns_Xi = (sum(d(1:(i-1)))+1):sum(d(1:i));
        JHk_i    = JHk(:, columns_Xi);
        JHkitA   = JHk_i'*A;
		alpha_i  = alphas(i);
		rho_i    = rhos(i);
		weight_i = d(i);  % weight when composited concave penalty is used
		lambda_i = lambdas(i)*weight_i*alpha_i/rho_i;
        
		for r=1:R,
            % form weights of the penalty according to previous sigma0_ir
		    sigma0_ir = Sigmas0(i,r);
            B_ir0     = B0(columns_Xi,r);
            omega_ir_g  = hfun_sg(sigma0_ir,gamma,1); % group level gradient
            omega_ir    = omega_ir_g*hfun_sg(abs(B_ir0),gamma,1); % individual level gradient
            
            % proximal operator of weighted L1 norm
            JHkitA_r  = JHkitA(:,r);
			lambda_ir = lambda_i*omega_ir;
			B_ir = sign(JHkitA_r) .* max(0, abs(JHkitA_r) - lambda_ir);
			B(columns_Xi,r) = B_ir;
			sigma_ir  = sum(hfun(abs(B_ir),gamma,1)); % sigma_{lr} = \sum g(blrj);
			Sigmas(i,r) = sigma_ir;
		end
    end
    
    % diagnostics
	Theta = ones(n,1)*mu + A*B';
	
	f_obj = 0;
	g_obj = 0;
    for i = 1:nDataSets,
        columns_Xi = (sum(d(1:(i-1)))+1):sum(d(1:i));
	    X_i        = X(:, columns_Xi);
	    W_i        = W(:, columns_Xi);
	    Theta_i    = Theta(:,columns_Xi);
	    alpha_i    = alphas(i);
		weight_i   = d(i);  % weight when composited concave penalty is used
		lambda_i   = lambdas(i)*weight_i;
    
	    log_partition = str2func([dataTypes{i} '_b']);
		W_it = W_i';
		Theta_it = Theta_i';
	    f_obj = f_obj + (1/alpha_i)*...
	        (fastTrace(W_it,log_partition(Theta_i)) - fastTrace(Theta_it,X_i));
	
	    for r = 1:R,
			sigma_ir    = Sigmas(i,r);
		    g_obj = g_obj + lambda_i*hfun(sigma_ir,gamma,1); % sigma_{lr} = \sum g(blrj)
        end
    end
	obj = f_obj + g_obj;
    
    % reporting
    out.hist_obj(k+1) = obj;   % objective + penalty
    out.f_obj(k+1)    = f_obj; % objective
    out.g_obj(k+1)    = g_obj; % penalty

	out.rel_obj(k)    = abs(obj0-obj)/abs(obj0); % relative change of loss function
	out.rel_Theta(k)  = norm(Theta0-Theta,'fro')^2/norm(Theta0,'fro')^2; % relative change of parameters
    
    % stopping checks
    if (out.rel_obj(k) < tol_obj); break; end;
	
	% remove the all zeros columns to simplify the computation and save memory
    if (threPath == 0),
        nonZeros_index = not(all(B==0,1));
        A = A(:,nonZeros_index);
        B = B(:,nonZeros_index);
	    Sigmas = Sigmas(:,nonZeros_index);
        [~,R]  = size(B);
    end
    
    % save previous results
    B0 = B; 
    Theta0  = Theta;
    Sigmas0 = Sigmas;	
    obj0    = obj;
end

 % variance explained for each data types
 varExpTotals = nan(1,nDataSets + 1);  % +1 is for the full data set
 varExpPCs    = nan(nDataSets + 1, R); % +1 is for the full data set
 X_full = nan(size(X)); % combine the quantitative data and the pseudo data of nonGaussian data
 for i = 1:nDataSets,
     columns_Xi = (sum(d(1:(i-1)))+1):sum(d(1:i));
     Xi    = X(:,columns_Xi);
	 mu_Xi = mu(1,columns_Xi);
     if not(strcmp(dataTypes{i}, 'Gaussian')),
         Xi = JHk(:,columns_Xi) + ones(n,1)*mu_Xi;
     end
	 X_full(:,columns_Xi) = Xi;
     
	 B_Xi  = B(columns_Xi,:);
     W_Xi  = W(:,columns_Xi);
	 
	 [tmp_total,tmp_PCs] = varExp_Gaussian(Xi,mu_Xi',A,B_Xi,W_Xi);
	 varExpTotals(1,i)   = tmp_total;
	 varExpPCs(i,:)      = tmp_PCs;
 end
 
 % variance explained ratio of the full data set
 weighted_X  = nan(size(X));
 weighted_mu = nan(size(mu));
 weighted_B  = nan(size(B));
 for i = 1:nDataSets,
     columns_Xi = (sum(d(1:(i-1)))+1):sum(d(1:i));
     Xi    = X(:,columns_Xi);
     mu_Xi   = mu(1,columns_Xi);
     if not(strcmp(dataTypes{i}, 'Gaussian')),
         Xi = JHk(:,columns_Xi) + ones(n,1)*mu_Xi;
     end
	 B_Xi    = B(columns_Xi,:);
	 alpha_i = alphas(i);
	 
	 weighted_X(:,columns_Xi)  = (1/sqrt(alpha_i))*Xi;
	 weighted_mu(1,columns_Xi) = (1/sqrt(alpha_i))*mu_Xi;
	 weighted_B(columns_Xi,:)  = (1/sqrt(alpha_i))*B_Xi;
 end
 
 [tmp_total,tmp_PCs] = varExp_Gaussian(X_full,mu',A,B,W);
 varExpTotals(1,end) = tmp_total;
 varExpPCs(end,:)    = tmp_PCs;
 
 % extract the strcuture index
 S = zeros(nDataSets, R);
 S(varExpPCs(1:nDataSets,:) > 0) = 1;

 out.iter = k;
 mu = mu';
 out.varExpTotals = varExpTotals;
 out.varExpPCs    = varExpPCs;
 out.Sigmas       = Sigmas;

end



