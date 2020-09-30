function model = truncatedEigenRegression(X,y,isXtX,opt)
%
% truncatedEigenRegression
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Implement linear regression y = X * B with a truncated eigenvalue
% regularisation of X' * X.
%
% Input arguments:
%
% X: data matrix of size [nPnts,p] or if isXtX == true, covariance matrix 
% X'*X) of size [p,p].
%
% y: data matrix of size [nPnts,q] or if isXtX == true, cross-covariance
% matrix (X'*y) of size [p,q].
%
% isXtX: boolean: whether X & y are observation (false) or covariance
% matrices (true).
%
% opt: structure containing regularisation options. In particular:
%
%   opt.nKeepDims: vector of integers indicating the number of dimensions
%   to keep (one model fitted per element).
%
%   opt.removeEig: options passed to 'LM_regEigen' for the initial
%   eigenvalue decomposition
%
% Output:
%
% model: structure containing:
%
% model.coeffs: matrix of size [p,q,numel(opt.nKeepDims)]: fitted
% coefficients
%
% model.nKeepDims: vectore of integers: number of dimensions of X'*X kept
% for each model
%
% model.eigenVals: sorted eigenvalues of X'*X
%
% See also: LM.fitLinearModel, LM.regEigen, LM.ridge
% ------
[~,p] = size(X);
[~,q] = size(y);

if ~isfield(opt,'removeEig')
    % use default tolerance in eigenvalue decomposition
    % (throw away degenerate dimensions)
    opt.removeEig = [];
end

% eigenvalue decomposition of X' * X & dropping numerical zeros
[V,S] = LM.regEigen(X, isXtX, opt.removeEig);

if isXtX
    % in this case y is actually X'*y precomputed!
    z = V' * y;
else
    z = V' * (X' * y); % KEEP the parenthesis!
end

nKeepDims = opt.nKeepDims;
nKeepDims = nKeepDims( nKeepDims <= numel(S) );
nOut = numel(nKeepDims);

coeffs = zeros(p,q,nOut);

for iOut = 1:nOut
    coeffs(:,:,iOut) = V(:,1:nKeepDims(iOut)) * ( z(1:nKeepDims(iOut),:) ./ S(1:nKeepDims(iOut)) );
end

% storing results
model = struct();
model.eigenVals = S;
model.nKeepDims = nKeepDims;
model.coeffs = coeffs;

end
%
%