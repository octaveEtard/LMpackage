function model = LM_fitLinearModel(X,y,opt)
%
% LM_fitLinearModel
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Implement linear regression y = X * B with regularisation of X' X.
%
% Input arguments:
%
% X: data matrix of size [nPnts,p] or if opt.accumulate == true, covariance
% matrix (X'*X) of size [p,p].
%
% y: data matrix of size [nPnts,q] or if opt.accumulate == true,
% cross-covariance matrix (X'*y) of size [p,q].
%
% opt: structure containing fitting options. In particular:
%
%   opt.accumulate: boolean: whether X & y are observation (false) or 
% covariance matrices (true).
%
%   opt.method: argument for the fitting algorithm. The algorithm used is
% defined by opt.method.name:
%
%       'trunc-eig-XtX': truncated eigenvalue regularisation of X' * X
%       'ridge-eig-XtX': ridge regularisation of X' * X
%
% See also: LM_ridge, LM_truncatedEigenRegression
%
% ------ sanity checks
assert( ~isempty(y) && ~isempty(X) );

switch opt.method.name
    
    case 'trunc-eig-XtX'
        model = LM_truncatedEigenRegression(X,y,opt.accumulate,opt.method);

    case 'ridge-eig-XtX'
        model = LM_ridge(X,y,opt.accumulate,opt.method);
end

end
%
%