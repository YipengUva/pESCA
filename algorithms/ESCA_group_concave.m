function [mu,A,B,S,out] = ESCA_group_concave(dataSets,dataTypes,lambdas,fun,opts)
% Expoential family simultaneous component analysis (ESCA) model with a group concave penalty
% on the loading matrix to induce structured sparsity on the loading matrix to seperate 
% the global, local common and distinct variation in multiple mixed types data sets.
%
% A Majorization-Minimization (MM) algorithm is used to fit the model. Missing value is allowed.
%
%   minimize - \sum_{l}^{L} (logP(X_{l}|Theta_{l}) + lambda_l*weight_l*(\sum_r g(type(B_{lr}))); 
%   s.t. Theta_l = 1*mu_l' + AB_l', l = 1...L; 
%            A'A = I;
%            1'A = 0;
%            type(B_{lr}) = ||B_{lr}||_2 % L2 norm;
%            type(B_{lr}) = ||B_{lr}||_1 % L1 norm;
%            type(B_{lr}) = \sum_{j} g(b_{lrj}) % composite concave penalty.
%
% Input:
%      dataSets: a cell array to hold the L data sets;
%      dataTypes: a cell array to hold the data types of the L data sets;
%      lambdas: a vector of values of the tuning parameter lambda;
%      fun: the name of the used concave penalty;
%      opts.
%           tol_obj: tolerance for relative change of hist_obj, default:1e-6;
%           maxit: max number of iterations, default: 1000;
%           gamma: hyper-parameter of penalty g(), default: 1;
%           R: the number of PCs, default: 50;
%           randomstart: initilization method, random (1), SCA(0), default: 0;
%           alphas: dispersion parameters of exponential dispersion families, default: 1. 
%           threPath: the option to generate thresholding path, default: 0;
%           type: 'L2', concave L2-norm penalty, only generate group sparsity;
%                 'L1', concave L1-norm penlaty, generate group sparsity, and slight elementwise sparsity;
%                 'biLevel', composite concave penalty, generate group and elementwise sparsity;
%           quiet: quiet==1, not show the progress when running the algorithm, default: 0.
%
% Output:
%       mu: offset term; 
%       A: score matrix;
%       B: loading matrix;
%       S: group sparse pattern of B;
%       out. 
%           iter: number of iterations;
%           hist_obj: joint loss function at each iteration;
%           f_obj: data fitting term at each iteration;
%           g_obj: regularization ata each iteration;
%           rel_obj: relative change of hist_obj at each iteration;
%           rel_Theta: relative change of Theta at each iteration;
%           varExpTotals: total variation explained of each data set;
%           varExpPCs: variation explained of each data set and each component;
%           Sigmas: the group length (the definition depends on the used input type).

% set the default type to be concave 2-norm penalty
if isfield(opts, 'type'), type = opts.type;  else type = 'L2'; end

% select proper algorithm according to the type
if strcmp(type, 'L2'),
    [mu,A,B,S,out] = ESCA_group_concave_L2(dataSets,dataTypes,lambdas,fun,opts);
elseif strcmp(type, 'biLevel'),
    [mu,A,B,S,out] = ESCA_group_concave_composite(dataSets,dataTypes,lambdas,fun,opts);
elseif strcmp(type, 'L1'),
    [mu,A,B,S,out] = ESCA_group_concave_L1(dataSets,dataTypes,lambdas,fun,opts);
end

end



