function [V,S] = regEigen(X,isXtX,reg)
%
% regEigen
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Implement eigenvalue decomposition : X' * X = V * S * V'
% and return V and S (as a vector) sorted by eigenval and with some
% regularisation.
%
% Input X will not be normalised in this function.
%
% Input arguments:
%
% X: data matrix of size [nPnts,p] or if isXtX == true, covariance matrix
% (X'*X) of size [p,p].
%
% isXtX: boolean indicating whether X is an observation (isXtX == false) or 
% covariance matrix (isXtX == true). Default is false.
%
% reg: structure containing regularisation options. Valid options are:
%
% reg.type =
%   'nKeep': keep reg.n dimensions
%	'var': keep dimensions accounting up to at most a fraction
%   reg.var of the total variance
% 	'tol': keep dimensions associated with an eigenvalue > reg.tol
%   'cond': keep dimensions up to a condition number < reg.cond
%   'ridge': transform eigenvalues by S -> S^2 / (S + reg.lambda * meanEigenvalue)
%
% Degenerate dimensions associated with null eigenvalues based on default
% tolerance will always be removed.
% ------
if nargin < 2
    isXtX = false;
end

if isXtX
    % input matrix is X' * X precomputed
    [V,S] = eig(X,'vector');
else
    [V,S] = eig(X'*X,'vector');
end

% eig does not sort the eigenvalues as svd !
[S,iSort] = sort(S,'descend');

% we'll remove at least these (numerical zeros)
rx = sum( eps(S(1)) < S);

% various criterion to choose how many components to keep
if 2 < nargin && ~isempty(reg)
    
    switch reg.type
        
        case 'nKeep'
            % keeping a number of dimensions
            rx = min( rx, reg.n );
            
        case 'var'
            % keeping dimensions accounting for a given fraction of variance
            rx = min( rx, sum( cumsum(S)./ sum(S) <= reg.var ) );
            
        case 'tol'
            % custom tolerance
            rx = min(rx, sum(reg.tol < S) );
            
        case 'cond'
            % condition number
            rx = min( rx, sum( (S(1)./S) < reg.cond ) );
            
        case 'ridge'
            % ridge
            S = S.^2 ./ ( S + reg.lambda * mean(S(1:rx)) );
            
    end
end

% rx = number of dimensions kept (rank of truncated X'*X)
% dropping dimensions & sorting V as S
V = V(:,iSort(1:rx));
S = S(1:rx);

end
%
%