function model = LM_ridge(X,y,isXtX,opt)
%
% LM_ridge
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Implement linear regression y = X * B with ridge regularisation of X' * X
%
% The regression is implemented through an eigenvalue decomposition
% of X' * X allowing to fit the model for a large number of regularisation
% parameters efficiently.
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
%   opt.lambda: vector of regularisation values to use
%
%   opt.normaliseLambda: boolean ; if true, the regularisation values
%   lambdas will be changed to lambda * (mean eigenvalues of X' * X)
%
%   opt.removeEig: options passed to 'LM_regEigen' for the initial
%   eigenvalue decomposition. If empty, use default tolerance.
%
% Output:
%
% model: structure containing:
%
% model.coeffs: matrix of size [p,q,numel(opt.lambda)]: fitted
% coefficients
%
% model.eigenVals: sorted eigenvalues of X'*X
%
% if opt.normaliseLambda: opt.meanEigenval : mean eigenvalue of X'*X that
% was used to normalise the regularisation coefficients.
%
% See also: LM_fitLinearModel, LM_regEigen, LM_truncatedEigenRegression
% ------

% !!! Assuming N > p here !
[~,p] = size(X);
[~,q] = size(y);

if ~isfield(opt,'removeEig')
    % use default tolerance in eigenvalue decomposition
    % (throw away degenerate dimensions)
    opt.removeEig = [];
end

% eigenvalue decomposition of X' * X & dropping numerical zeros
[V,S] = LM_regEigen(X, isXtX, opt.removeEig);


%% ------ ridge ------
lambda = opt.lambda;

if opt.normaliseLambda
    % normalising lambdas
    nL = mean(S);
    lambda = nL * lambda;
end

if isXtX
    % in this case y is actually X'*y precomputed!
    z = V' * y;
else
    z = V' * (X' * y); % KEEP the parenthesis!
end

nLambda = length(lambda);
coeffs = zeros(p,q,nLambda);

for iL = 1:nLambda
    coeffs(:,:,iL) = V * ( z ./ (S + lambda(iL)) );
end

% storing results
model = struct();

model.eigenVals = S;
model.coeffs = coeffs;

if opt.normaliseLambda
    model.meanEigenval = nL;
end

end
%
%