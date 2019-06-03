%% Demo of the P-ESCA model on Gaussian-Gaussian-Gaussian data sets
% This doc is going to show how to simulate Gaussian-Gaussian-Gaussian
% (G-G-G) data sets with underlying global, local common and distinct 
% structures according to the ESCA model. After that, these data sets are
% used to illustrate how to construct a P-ESCA model and the corresponding model
% selection process. In this doc, we will use the simulated parameters to
% evluate the model selection process. The documents of the used functions 
% can be found by using the command 'doc function name'.

%% Add current folder to the path
clear all
current_fold = pwd;
addpath(genpath(current_fold));

%% How to simulate proper G-G-G data sets
% G-G-G data sets are simulated according to the ESCA model. The number of
% samples is set to 100; the number of the variables in the
% three data sets are 1000, 500 and 100; the SNRs in simulating the global,
% local common and distinct structures are set to 1; the noise levels (the
% square root of the dispersion parameter $\alpha$) are set to 1; the
% sparse ratio, which is the proportion of 0s in the simulated loading
% matrix, is set to 0; ntimes_noise, which indicates that the singular values 
% of any specific structure are at least ntimes larger than the singular 
% value of the corresponding noise term, is set to 0.
% After that, the signals, which are the singular values of the simulated structures,
% and the noise, which are the singular values of the corresponding noise terms, are
% shown in Fig.~1. The variation explained ratios of each component for every
% single data set are shown in Fig.~2. 

%
% parameters used in the simulation 
n  = 100; % samples
ds = [1000, 500, 100]; % number of variables of each data set
noises = [1,1,1]; % noise levels
sparse_ratio = 0; % sparse ratio of the loading matrix
ntimes_noise = 0; % require signal strength is ntimes larger than the noise

SNRgc = 1;         % SNR of the global common structure
SNRlc = [1, 1, 1]; % SNRs of the local common structures
SNRd  = [1, 1, 1]; % SNRs of the distinct structures

% data simulation process
seed = 1234; % set a seed to reproduce the example, can be not exist
[dataSimulation] = ...
     dataSimulation_GGG(n,ds,SNRgc,SNRlc,SNRd,...
     noises,ntimes_noise,sparse_ratio,seed);

X1 = dataSimulation.X(:,1:ds(1));
X2 = dataSimulation.X(:,(ds(1)+1):(sum(ds(1:2))));
X3 = dataSimulation.X(:,(sum(ds(1:2))+1):end);

% characterize the singular values of the signal and noise
subTitles ={'C123','C12','C13','C23','D1','D2','D3'};

% form the block sparse pattern
blocks_sparse_index    = cell(7,1);
blocks_sparse_index{1} = 1:sum(ds);        % C123
blocks_sparse_index{2} = 1:sum(ds(1:2));   % C12
blocks_sparse_index{3} = [1:ds(1), (sum(ds(1:2))+1):sum(ds)]; % C13
blocks_sparse_index{4} = (sum(ds(1))+1):sum(ds); % C23
blocks_sparse_index{5} = 1:sum(ds(1)); % D1
blocks_sparse_index{6} = (ds(1)+1):(sum(ds(1:2))); % D2
blocks_sparse_index{7} = (sum(ds(1:2))+1):sum(ds); % D3

figure;
for i = 1:7
    index_factors = (3*(i - 1)+1):3*i;
    Theta_factors = dataSimulation.U_simu(:,index_factors)*...
                    diag(dataSimulation.D_simu(index_factors,1))*...
                    dataSimulation.V_simu(:,index_factors)';
    index_variables = blocks_sparse_index{i};
    Theta_factors = Theta_factors(:,index_variables);
    E_factors     = dataSimulation.E_simu(:,index_variables);
    
    subplot(1,7,i);
    [~,signal_factors,~] = fastSVD(Theta_factors,3);
    signal_factors = diag(signal_factors);
    [~,noise_factors,~]  = fastSVD(E_factors,3);
    noise_factors = diag(noise_factors);
    plot(signal_factors, '-o'); hold on;
    plot(noise_factors, '-o'); hold on;
    title(subTitles{i});
	if (i==1), ylabel('singular value'); xlabel('components'); end;
end
snapnow
disp(['Fig.~1: The singular values of the signal (blue dots),'...
                                       'of the noise (red dots).']);

% characterize the variation explained by each component for every data set
figure
elementLabels = arrayfun(@(x){sprintf('%0.1f%%',x)},...
                              dataSimulation.varExpPCs_simu(1:3,:));
colLabels     = arrayfun(@(x){sprintf('%0.2f%%',x)},...
                              dataSimulation.varExpTotals_simu);
heatmap(dataSimulation.varExpPCs_simu(1:3,:), ...
        1:size(dataSimulation.varExpPCs_simu,2), colLabels,...
        elementLabels, 'ShowAllTicks', true, 'Colormap', 'red',...
        'Colorbar', true);
title('true variation explained for each data set');
snapnow
disp('Fig.~2: Variation explained ratios for the three simulated data sets.'); 

%% How to estimate the dispersion parameter $\alpha$ for each data set
% The dispersion parameters are estimated using the PCA model.
% Details can be found in the supplementary information of the paper. The
% CV error plots of the the rank selection in the PCA model of
% these three data sets are shown in Fig.~3, Fig.~4, Fig.~5 respectively.
% The estimated dispersion parameters for the three data sets are shown in 
% alphas.

% 
% data sets and data types
dataSets  = {X1,X2,X3};
dataTypes = {'Gaussian', 'Gaussian', 'Gaussian'}; 

Rs = 5:15; % searching range
K  = 3;    % K times repetations
alphas = zeros(1,3); % dispersion parameters
ranks_estimation = zeros(3,K);

for (i=1:3),
    X_i = dataSets{i};
    dataType_i = dataTypes{i};
    if strcmp(dataType_i, 'Gaussian'),
        [est_mean,~,R_CV] = alpha_estimation(X_i,K,Rs);
        alphas(i) = est_mean;
        ranks_estimation(i,:) = R_CV;
    else
        alphas(i) = 1;
    end
    
    snapnow
    disp(sprintf(['Fig.~' num2str(i+2) ': CV error based model selection'...
                                  ' for the ' num2str(i) '-th data set']));
end
alphas
% when estimated alphas have large difference or very small, it is better to 
% scale the data set by sqrt(alpha) and then set alphas to be 1.

%% How to set the parameters for the P-ESCA model
% The following parameters are used for the P-ESCA model selection. The
% meaning of these parameters can be found in the document of the function.

fun = 'GDP'; gamma = 1; % penalty function
% fun = 'lp'; gamma = 0.5; % Lq penalty, q = 0.5
% fun = 'lp'; gamma = 1;   % Lq penalty, q = 1, group lasso penalty

opts = [];
opts.gamma = gamma;  % hyper-parameter for the used penalty
opts.random_start = 0;
opts.tol_obj = 1e-6; % stopping criteria 
opts.maxit   = 500;
opts.alphas  = alphas;
opts.R    = 50;    % components used 
opts.threPath = 1; % generaint thresholding path or not
opts.type = 'L2';  % induce group sparsity
%opts.type = 'L1'; % induce group sparsity + slight elementwise sparsity
%opts.type = 'biLevel'; % induce group sparsity + elementwise sparsity

%opts.quiet = 1; % don't show the progress when running the algorithm

%% The model selection process of the P-ESCA model
% At first, we need to find a proper searching range of $\lambda$ for the 
% model selection. After that, 30 values of $\lambda$ are selected from this
% searching range equally in log-space. Then, the developed CV error based
% model selection procedure is used to find a proper model. Fig.~6 shows
% how regularization strength $\lambda$ effects the CV errors (top left),
% the the RMSEs (top right), the RV coefficients in estimating the common
% structures (bottom left), and distinct structures (bottom right) when the
% P-ESCA model with a group GDP penalty is used. Fig.~7
% shows how $\lambda$ effects the rank estimation of the common and
% distinct structures. The red cross marker indicates the point corresponding 
% to the minimum CV error. Fig.~8 shows the thresholding paths during the
% model selection process. The red line indicates the selected cutting
% point.

% find a proper searching range of lambda
opts.maxit  = 50;
lambda_test = 20;
lambdas     = lambda_test*ones(1,3);

[~,~,~,~,out] = ESCA_group_concave(dataSets,dataTypes,lambdas,fun,opts);
out.varExpPCs(1:3,:)

%% specify the searching range of lambda
opts.maxit = 500;
K = 1;
n_tries  = 30;
lambda_i = 1;
lambda_t = 500;
lambdas_md = logspace(log10(lambda_i),log10(lambda_t),n_tries);

% model selection process
[cvErrors_mat,ranks_mat,RVs_mat,RMSEs_mat,inits,outs] = ...
                       ESCA_modelSelection_KCV_fullInfo(ds,dataTypes,...
                       K,lambdas_md,dataSimulation,fun,opts);
                    
%% how CV errors, RMSEs and RVs change with respect to the value of lambda
% index out the point corresponding to the minimum CV error
[~, index_min] = min(cvErrors_mat(:,1));
figure
subplot(2,2,1); % CV errors
plot(log10(lambdas_md), cvErrors_mat,'-o'); hold on;
plot(log10(lambdas_md(index_min)), cvErrors_mat(index_min,:),...
                      'rx','MarkerSize',10); hold on;
xlabel('log_{10}(\lambda)');
ylabel('CV error');
title('CV errors');
xlim([min(log10(lambdas_md)),max(log10(lambdas_md))]);
legend('X','X_1','X_2','X_3','Location','best')
subplot(2,2,2); % RMSEs
plot(log10(lambdas_md),RMSEs_mat,'-o'); hold on;
plot(log10(lambdas_md(index_min)), RMSEs_mat(index_min,:),...
            'rx','MarkerSize',10); hold on;
xlabel('log_{10}(\lambda)');
ylabel('RMSE');
title('RMSEs');
xlim([min(log10(lambdas_md)),max(log10(lambdas_md))]);
legend('\Theta','\Theta_1','\Theta_2','\Theta_3','\mu','Location','best')
subplot(2,2,3); % RV ceofficients for the common structures
plot(log10(lambdas_md),RVs_mat(:,1:4),'-o'); hold on;
plot(log10(lambdas_md(index_min)), RVs_mat(index_min,1:4),...
                'rx','MarkerSize',10); hold on;
xlabel('log_{10}(\lambda)');
ylabel('RV coefficients');
xlim([min(log10(lambdas_md)),max(log10(lambdas_md))]);
legend('C123','C12','C13','C23','Location','best')
title('common structures');
subplot(2,2,4); % RV ceofficients for the distinct structures
plot(log10(lambdas_md),RVs_mat(:,5:end),'-o'); hold on;
plot(log10(lambdas_md(index_min)), RVs_mat(index_min,5:end),...
                             'rx','MarkerSize',10); hold on;
xlabel('log_{10}(\lambda)');
ylabel('RV coefficients');
xlim([min(log10(lambdas_md)),max(log10(lambdas_md))]);
legend('D1','D2','D3','Location','best')
title('distinct structures');

snapnow
disp(sprintf(['Fig.~6: How CV errors, RMSEs and RV coefficients \n ' ...
    'change during the model selection process of the P-ESCA model.']));

% how the rank estimation of any specific structure changes with respect
% to the value of $\lambda$.
figure
subplot(1,2,1)
plot(log10(lambdas_md),ranks_mat(:,1:4),'-o'); hold on; 
plot(log10(lambdas_md(index_min)), ranks_mat(index_min,1:4),...
                             'rx','MarkerSize',10); hold on;
xlabel('log_{10}(\lambda)');
ylabel('estimated ranks');
title('estimated ranks');
legend('C123','C12','C13','C23','Location','best')
subplot(1,2,2)
plot(log10(lambdas_md),ranks_mat(:,5:end),'-o'); hold on;
plot(log10(lambdas_md(index_min)), ranks_mat(index_min,5:end),...
                               'rx','MarkerSize',10); hold on;
xlabel('log_{10}(\lambda)');
ylabel('estimated ranks');
title('estimated ranks');
legend('D1','D2','D3','Location','best')
snapnow
disp(sprintf(['Fig.~7: How rank estimations of the estimated structures \n' ... 
    'change during the model selection process of the P-ESCA model.']));

% thresholding pathes
% index out the Sigmas to generate thresholding path
Sigmas_threPath = zeros(3,opts.R,n_tries);
if strcmp(opts.type, 'L2'),
    weights = sqrt(ds');
else
    weights = ds';
end

% scaling the Sigmas according to the weights
for (i=1:n_tries),
    out_tmp = outs{i,1};
    Sigmas_tmp = out_tmp.Sigmas;
    Sigmas_threPath(:,:,i) = Sigmas_tmp ./ (sqrt(ds') *ones(1,50));
end
maxLen = max(Sigmas_threPath(:))*1.1;

% thresholding profiles plot
figure
subplot(1,3,1)
lambdas_md_X1 = lambdas_md;
plot(log10(lambdas_md_X1), (squeeze(Sigmas_threPath(1,:,:)))','-'); hold on;
plot([log10(lambdas_md_X1(index_min)),log10(lambdas_md_X1(index_min))],...
     [0,  maxLen],'r');
title('X1');
ylim([0, maxLen]);
xlim([min(log10(lambdas_md_X1)), max(log10(lambdas_md_X1))]);
xlabel('log_{10}(\lambda)');
ylabel('scaled L_2 norm');
subplot(1,3,2)
lambdas_md_X2 = lambdas_md;
plot(log10(lambdas_md_X2), (squeeze(Sigmas_threPath(2,:,:)))','-'); hold on;
plot([log10(lambdas_md_X2(index_min)),log10(lambdas_md_X2(index_min))],...
     [0, maxLen],'r');
title('X2');
xlim([min(log10(lambdas_md_X2)), max(log10(lambdas_md_X2))]);
xlabel('log_{10}(\lambda)');
ylim([0, maxLen]);
subplot(1,3,3)
lambdas_md_X3 = lambdas_md;
plot(log10(lambdas_md_X3), (squeeze(Sigmas_threPath(3,:,:)))','-'); hold on;
plot([log10(lambdas_md_X3(index_min)),log10(lambdas_md_X3(index_min))],...
     [0, maxLen],'r');
title('X3');
ylim([0, maxLen]);
xlim([min(log10(lambdas_md_X3)), max(log10(lambdas_md_X3))]);
xlabel('log_{10}(\lambda)');
snapnow
disp('Fig.~8: Thresholding pathes of the model selection of the P-ESCA model.');

%% How to fit the final model
% After selecting the model with the minimum CV error, the selected model is
% re-fitted on the full data set with the selected value of $\lambda$ and
% the outputs of the selected model as the initialization. The RV
% coefficients in estimating the global common (C123), local common (C12,
% C13, C23) and distinct structures (D1, D2, D3) are shown as RVs_structures. 
% And the corresponding rank estimations of these structures are shown as 
% Ranks_structures. The RMSEs in estimating the simulated
% $\mathbf{\Theta}$, $\mathbf{\Theta}_1$, $\mathbf{\Theta}_2$ and
% $\mathbf{\Theta}_3$ and $\mu$ are shown as RMSEs_parameters. Fig.~9 shows
% the variation explained ratios computed using the estimated parameters for 
% the three data sets.

%
% parameters used to fit the final model
opts_inner = inits{index_min,1};
opts_inner.tol_obj = 1e-8;
opts_inner.maxit   = 5000;
opts_inner.threPath = 0;

lambdas = lambdas_md(index_min)*ones(1,length(dataSets));
[mu,A,B,S,out] = ESCA_group_concave(dataSets,dataTypes,lambdas,fun,opts_inner);

%  evaluate the performance of the final model
[RVs_structures,Ranks_structures,RMSEs_parameters] = ...
    ESCA_evaluation_metrics(mu,A,B,S,ds,dataSimulation);
RVs_structures
Ranks_structures
RMSEs_parameters

%% estimated variation explainted ratios
% index out different structures
C123_index = sum(S,1)==3;
C12_index = sum(S([1,2],:),1)== 2 - C123_index;
C13_index = sum(S([1,3],:),1)== 2 - C123_index;
C23_index = sum(S([2,3],:),1)== 2 - C123_index;
D_index = sum(S(:,:),1) == 1;
D1_index = D_index & sum(S(1,:),1) == 1;
D2_index = D_index & sum(S(2,:),1) == 1;
D3_index = D_index & sum(S(3,:),1) == 1;

varExpPCs    = out.varExpPCs(1:3,:);
varExpTotals = out.varExpTotals(1:3);
varExpPCs_plots = [varExpPCs(:,C123_index), varExpPCs(:,C12_index),...
                   varExpPCs(:,C13_index),varExpPCs(:,C23_index), ...
                   varExpPCs(:,D1_index), varExpPCs(:,D2_index),...
                   varExpPCs(:,D3_index)];
               
% variation explained ratios
figure
elementLabels = arrayfun(@(x){sprintf('%0.1f%%',x)}, varExpPCs_plots);
colLabels     = arrayfun(@(x){sprintf('%0.2f%%',x)}, varExpTotals);
heatmap(varExpPCs_plots,  1:size(varExpPCs,2), colLabels, elementLabels,...
        'ShowAllTicks', true, 'Colormap', 'red',...
        'Colorbar', true);
title('estimated variation explained');
xlabel('components');
snapnow
disp('Fig.~9: The estimated variation explained ratios.');

