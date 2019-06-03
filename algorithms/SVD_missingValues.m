function [U,S,V,out] = SVD_missingValues(X,R,opts)
% SVD algorithm with the option for missing values
%   minimize 0.5*||W.*(X-Z)||_F^2; 
%   s.t. Z = U*S*V';
%   W: weighting matrix for missing values
%
% Input:
%      X: m*n matrix
%      opts.
%           tol_obj: tolerance for relative change of function value, default:1e-4
%           maxit: max number of iterations, default: 1000
%
% Output:
%       U, S, V: same as standard svd algorithm  
%       out. 
%           iter: number of iterations
%           hist_obj: objective value at each iteration
%           relerr_obj: relative change of objective function at each iteration

% parameters and defaults
if isfield(opts, 'tol_obj'),  tol_obj = opts.tol_obj;  else tol_obj = 1e-4;  end
if isfield(opts, 'maxit'),    maxit   = opts.maxit;    else maxit   = 1000;  end

% form weighting matrix
Wc = isnan(X);
W  = 1 - Wc;      % weighting matrix
X(isnan(X)) = 0;  % remove missing elements

% initialization
% initial parameters
if(isfield(opts, 'U0'))  % using inputted initialization
    U0 = opts.U0;
    S0 = opts.S0;
	V0 = opts.V0;
    Z0 = U0(:,1:R)* S0(1:R,1:R) * V0(:,1:R)';
else 
    [U0,S0,V0] = fastSVD(X,R);
    Z0 = U0*S0*V0';	
end

% initial value of loss function
obj0 = 0.5*norm(W.*(X-Z0),'fro')^2;
out.hist_obj(1) = obj0;

% iterations
for k = 1:maxit
    %fprintf('\n%d th iteration',k);
	
	% form Xtilde 
	Xtilde = W.*X + Wc.*Z0;
	
	% update Z
	[U,S,V] = fastSVD(Xtilde, R);
	Z = U*S*V';
	
	% new objective value
	obj = 0.5*norm(W.*(X-Z),'fro')^2;
	 
    % reporting
    out.hist_obj(k+1) = obj;
	out.relerr_obj(k) = (obj0-obj)/(obj0); 
	
    % stopping checks
    if ((out.relerr_obj(k)<tol_obj)); break; end;
    
    % save previous results
    Z0   = Z; 
    obj0 = obj; 
end
 
 out.iter = k; 
 %fprintf('Convergenced in %d iterations. \n',k);  % report # of iterations
end
    
