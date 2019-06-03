function [ out ] = obj_ls_gradient(X,Theta)
% gradient of least square loss
% X: data matrix
% Theta: parameter matrix 
out = Theta - X;
end

