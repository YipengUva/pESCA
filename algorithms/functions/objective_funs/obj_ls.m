function [ out ] = obj_ls(X,Theta,P)
% object function of least square loss
% X: data matrix
% Theta: parameter matrix 
out = 0.5*norm(P.*(X-Theta),'fro')^2;
end

