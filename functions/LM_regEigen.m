function [V,S] = LM_regEigen(X,isXtX,reg)
%
% LM_regEigen
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Implement eigenvalue decomposition : X' * X = V * S * V'
% and return V and S (as a vector) sorted by eigenval and with some
% regularisation.
%
% Input X will not be normalised in this function
%
% ------
if nargin < 2
    isXtX = false;
end

if isXtX
    % input matrix is X' * X p precomputed
    [V,S] = eig(X,'vector');
else
    [V,S] = eig(X'*X,'vector');
end

% eig does not sort the eigenvalues as svd !
[S,iSort] = sort(S,'descend');

% we'll remove at least these (numerical zeros)
rx = sum(S > eps(S(1)));

% various criterion to choose how many components to keep
if 2 < nargin && ~isempty(reg)
    
    switch reg.type
        
        case 'nKeep'
            % keeping a number of dimensions
            rx = min( rx, reg.n );
            
        case 'var'
            % keeping dimensions accounting for a given fraction of variance
            rx = min( rx, find( cumsum(S)./ sum(S) > reg.var,1)-1 );
            
        case 'tol'
            % custom tolerance
            rx = min(rx, sum(S > reg.tol) );
            
        case 'cond'
            % condition number
            rx = min( rx, find(S(1)./S > reg.cond,1)-1 );
            
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