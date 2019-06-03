function [dataSimulation] = dataSimulation_GBB_pro(n,ds,SNRgc,SNRlc,SNRd,noises,margProb,sparse_ratio)

% This function will be used to simulate three blocks
% Gaussian-Bernoulli_Bernoulli data sets

% parameters used in the simulation
sumd = sum(ds);   % sum of variables
Rgc  = 3;         % number of PCs global common 
Rlc  = [3, 3, 3]; % number of PCs local common
Rd   = [3, 3, 3]; % number of PCs distint 
sumR = Rgc + sum(Rlc) + sum(Rd);

% SNRs used in the data simulation process
SNRs = [SNRgc, SNRlc, SNRd];

% simulate the offset term
mu_1    = randn(ds(1),1);
mu_23   = logit(betarnd(round(n*margProb)+1, n - round(n*margProb)+1, sum(ds(2:3)),1));
mu_simu = [mu_1;mu_23];

% the simulation of U, D_pre, V_pre
U_pre = mvnrnd(zeros(1,sumR),eye(sumR),n);
[U_simu,~,~] = svds(U_pre - ones(n,1)*mean(U_pre,1),sumR);
D_pre        = abs(1+0.5*randn(sumR,1));
[V_pre,~]    = qr(mvnrnd(zeros(1,sumR),eye(sumR),sumd),0);

% modify V_pre to have the predifed structure
% the target group sparsity
S_simu = [repmat([1;1;1],1,Rgc),...
          repmat([1;1;0],1,Rlc(1)),...
          repmat([1;0;1],1,Rlc(2)),...
          repmat([0;1;1],1,Rlc(3)),...
          repmat([1;0;0],1,Rd(1)),...
          repmat([0;1;0],1,Rd(2)),...
          repmat([0;0;1],1,Rd(3))];
      
% modify V_pre to have the target group sparsity
V_simu = nan(sumd,sumR);
for i = 1:3
    columns_Xi = (sum(ds(1:(i-1)))+1):sum(ds(1:i));
    for r = 1:sumR
        V_simu(columns_Xi,r) = sign(S_simu(i,r))*V_pre(columns_Xi,r);
    end
end

% impose individual sparsity
for i = 1:3
    columns_Xi = (sum(ds(1:(i-1)))+1):sum(ds(1:i));
    for r = 1:sumR
        V_ir = V_simu(columns_Xi,r);
        
        % individual sparsity: randomly set sparse_ratio% elements to be 0
        index_tmp =  randsample(1:ds(i), ceil(sparse_ratio*ds(i)));
        V_ir(index_tmp,:) = 0;
        
%         % individual sparsity: set the smallest elements to be 0
%         abs_V_ir = abs(V_ir);
%         cut_tmp  = quantile(abs_V_ir, sparse_ratio);
%         V_ir(abs_V_ir < cut_tmp) = 0;
        
        V_simu(columns_Xi,r) = V_ir;        
    end
end

% noise follows standard normal distribution
SLogistic  = makedist('Logistic');
E_simu = [noises(1)*randn(n,ds(1)),...
          random(SLogistic,n,ds(2)),...
          random(SLogistic,n,ds(3))];
          
% Modify the D_pre as C.*D to satisfy the pre-defined SNR
C_factors = nan(sumR,1);

for i = 1:length(SNRs)
    index_factors = (3*(i - 1)+1):3*i;
    Theta_pre = U_simu(:,index_factors)*diag(D_pre(index_factors,1))*V_simu(:,index_factors)';
    index_variables = sum(abs(Theta_pre),1) > 0;
    Theta_pre = Theta_pre(:,index_variables);
    E_pre     = E_simu(:,index_variables);
    
    % compute the proper scale for SNR
    C_factors(index_factors,1) = sqrt(SNRs(1,i))*norm(E_pre,'fro')/norm(Theta_pre,'fro');
end

D_simu = C_factors.*D_pre;

% simulate Theta_simu
Theta_simu = ones(n,1)*mu_simu' + U_simu*(diag(D_simu))*V_simu';

% simulate Gaussian-Bernoulli-Bernoulli
X1  = Theta_simu(:,1:ds(1)) + E_simu(:,1:ds(1));

Theta23_simu = Theta_simu(:,(ds(1)+1):end);
X23 = sign(fai_logistic(Theta23_simu) - rand(n,sum(ds(2:3))));
X23 = (0.5).*(X23 + 1);
X  = [X1, X23];

% latent data sets
Xstar = Theta_simu + E_simu;

% variance explained for each data types
B_simu = V_simu*diag(D_simu);
varExpTotals_simu = nan(1,4);
varExpPCs_simu    = nan(4, sumR);
X_full = nan(size(X));

for i = 1:3,
    columns_Xi = (sum(ds(1:(i-1)))+1):sum(ds(1:i));
    Xi    = Xstar(:,columns_Xi);
    mu_Xi = mu_simu(columns_Xi,1);
    B_Xi  = B_simu(columns_Xi,:);
    W_Xi  = ones(size(Xi));
	X_full(:,columns_Xi) = Xi;
    
	[tmp_total,tmp_PCs]    = varExp_Gaussian(Xi,mu_Xi,U_simu,B_Xi,W_Xi);
	varExpTotals_simu(1,i) = tmp_total;
	varExpPCs_simu(i,:)    = tmp_PCs;
end
 
% variance explained ratio combine all the data sets
  W = ones(size(X));
 [tmp_total,tmp_PCs] = varExp_Gaussian(X_full,mu_simu,U_simu,B_simu,W);
 varExpTotals_simu(1,end) = tmp_total;
 varExpPCs_simu(end,:)    = tmp_PCs;
 
 % save the results
 % save X, Theta, E
 dataSimulation.X = X;
 dataSimulation.Theta_simu = Theta_simu;
 dataSimulation.E_simu     = E_simu;
 dataSimulation.mu_simu    = mu_simu;
 
 % save mu, U, D, V
 dataSimulation.mu_simu = mu_simu;
 dataSimulation.U_simu  = U_simu;
 dataSimulation.D_simu  = D_simu;
 dataSimulation.V_simu  = V_simu;
 
 % save true variation explained
 dataSimulation.varExpTotals_simu = varExpTotals_simu;
 dataSimulation.varExpPCs_simu    = varExpPCs_simu;

end

