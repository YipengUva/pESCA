function [RVs_structures,Ranks_structures,RMSEs_parameters] = ...
    ESCA_evaluation_metrics(mu,A,B,S,ds,dataSimulation)

% index out the simulated structures
X          = dataSimulation.X;
Theta_simu = dataSimulation.Theta_simu;
mu_simu    = dataSimulation.mu_simu;
U_simu     = dataSimulation.U_simu;
D_simu     = dataSimulation.D_simu;
V_simu     = dataSimulation.V_simu;

[n,~] = size(X);

% simulated parameters
% Theta1, Theta2, Theta3
Theta1_simu = Theta_simu(:,1:ds(1));
Theta2_simu = Theta_simu(:,(ds(1)+1):(sum(ds(1:2))));
Theta3_simu = Theta_simu(:,(sum(ds(1:2))+1):end);

% C123
i = 1; index_factors = (3*(i - 1)+1):3*i;
C123_simu = U_simu(:,index_factors)*diag(D_simu(index_factors,1))*V_simu(:,index_factors)';
% C12
i = 2; index_factors = (3*(i - 1)+1):3*i;
C12_simu  = U_simu(:,index_factors)*diag(D_simu(index_factors,1))*V_simu(:,index_factors)';
C12_simu  = C12_simu(:,sum(abs(C12_simu),1) > 0);
% C13
i = 3; index_factors = (3*(i - 1)+1):3*i;
C13_simu  = U_simu(:,index_factors)*diag(D_simu(index_factors,1))*V_simu(:,index_factors)';
C13_simu  = C13_simu(:,sum(abs(C13_simu),1) > 0);
% C23
i = 4; index_factors = (3*(i - 1)+1):3*i;
C23_simu  = U_simu(:,index_factors)*diag(D_simu(index_factors,1))*V_simu(:,index_factors)';
C23_simu  = C23_simu(:,sum(abs(C23_simu),1) > 0);
% D1
i = 5; index_factors = (3*(i - 1)+1):3*i;
D1_simu  = U_simu(:,index_factors)*diag(D_simu(index_factors,1))*V_simu(:,index_factors)';
D1_simu  = D1_simu(:,sum(abs(D1_simu),1) > 0);
% D2
i = 6; index_factors = (3*(i - 1)+1):3*i;
D2_simu  = U_simu(:,index_factors)*diag(D_simu(index_factors,1))*V_simu(:,index_factors)';
D2_simu  = D2_simu(:,sum(abs(D2_simu),1) > 0);
% D3
i = 7; index_factors = (3*(i - 1)+1):3*i;
D3_simu  = U_simu(:,index_factors)*diag(D_simu(index_factors,1))*V_simu(:,index_factors)';
D3_simu  = D3_simu(:,sum(abs(D3_simu),1) > 0);

% model evulation
% using simulated parameters
Theta_Hat  = ones(n,1)*mu' + A*B';
Theta1_Hat = Theta_Hat(:,1:ds(1));
Theta2_Hat = Theta_Hat(:,(ds(1)+1):(sum(ds(1:2))));
Theta3_Hat = Theta_Hat(:,(sum(ds(1:2))+1):end);

C123_index = sum(S,1)==3;
C123_Hat = A(:,C123_index)*B(:,C123_index)';
        
C12_index = sum(S([1,2],:),1)== 2 - C123_index;
C12_Hat = A(:,C12_index)*B(:,C12_index)';
C12_Hat = C12_Hat(:,sum(abs(C12_Hat),1) > 0);

C13_index = sum(S([1,3],:),1)== 2 - C123_index;
C13_Hat = A(:,C13_index)*B(:,C13_index)';
C13_Hat = C13_Hat(:,sum(abs(C13_Hat),1) > 0);

C23_index = sum(S([2,3],:),1)== 2 - C123_index;
C23_Hat = A(:,C23_index)*B(:,C23_index)';
C23_Hat = C23_Hat(:,sum(abs(C23_Hat),1) > 0);

D_index = sum(S(:,:),1) == 1;
D1_index = D_index & sum(S(1,:),1) == 1;
D2_index = D_index & sum(S(2,:),1) == 1;
D3_index = D_index & sum(S(3,:),1) == 1;

D1_Hat = A(:,D1_index)*B(:,D1_index)';
D1_Hat = D1_Hat(:,sum(abs(D1_Hat),1) > 0);
D2_Hat = A(:,D2_index)*B(:,D2_index)';
D2_Hat = D2_Hat(:,sum(abs(D2_Hat),1) > 0);
D3_Hat = A(:,D3_index)*B(:,D3_index)';
D3_Hat = D3_Hat(:,sum(abs(D3_Hat),1) > 0);
        
% RV coefficients of estimated structures
RV_C123 = RV_modified_bda(C123_simu, C123_Hat);
RV_C12  = RV_modified_bda(C12_simu, C12_Hat);
RV_C13  = RV_modified_bda(C13_simu, C13_Hat);
RV_C23  = RV_modified_bda(C23_simu, C23_Hat);
RV_D1   = RV_modified_bda(D1_Hat, D1_simu);
RV_D2   = RV_modified_bda(D2_Hat, D2_simu);
RV_D3   = RV_modified_bda(D3_Hat, D3_simu);
RVs_structures = [RV_C123,RV_C12,RV_C13,RV_C23,RV_D1,RV_D2,RV_D3];

% ranks of estimated structures
Ranks_structures = [sum(C123_index),sum(C12_index),sum(C13_index),...
              sum(C23_index),sum(D1_index),sum(D2_index),sum(D3_index)];
                      
% RV coefficients of estimated Theta
RMSE_Theta1 = norm(Theta1_Hat-Theta1_simu,'fro')^2/norm(Theta1_simu,'fro')^2;
RMSE_Theta2 = norm(Theta2_Hat-Theta2_simu,'fro')^2/norm(Theta2_simu,'fro')^2;
RMSE_Theta3 = norm(Theta3_Hat-Theta3_simu,'fro')^2/norm(Theta3_simu,'fro')^2;
RMSE_Theta  = norm(Theta_Hat-Theta_simu,'fro')^2/norm(Theta_simu,'fro')^2;
RMSE_mu     = norm(mu-mu_simu,'fro')^2/norm(mu_simu,'fro')^2;
RMSEs_parameters = [RMSE_Theta, RMSE_Theta1, RMSE_Theta2, RMSE_Theta3, RMSE_mu];

end

