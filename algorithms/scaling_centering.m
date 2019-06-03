function [X_new] = scaling_centering(X,scale,center)
% function to do column centering and scaling of a matrix
% center = 1 or 0, centering the X or not
% scale   = 1 or 0, scaling the X by the emperical std

% default parameters
% only scaling X
if (nargin<3), center = 0; end;
if (nargin<2), scale  = 1; end;

% number of samples
[n,~]  = size(X);

% centering X or not
if (center == 1)
    X_mean = mean(X,1,'omitnan');
    X_tmp  = X - ones(n,1)*X_mean;
else
    X_tmp  = X;
end

% scaling X or not
if (scale == 1)
    X_std = std(X,1,'omitnan');
    X_new = X_tmp./(ones(n,1)*X_std);
else
    X_new = X_tmp;
end

end

