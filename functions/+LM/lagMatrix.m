function X = lagMatrix(x, nLags)
%
% lagMatrix
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Given a matrix x with [nPnts,nCol] = size(x), and nLags integer,
% n = nPnts - nLags + 1, make matrix X of size [n,nLags * nCol]:
%
% X(:,(1:nLags) + (iCol-1) * nLags) =
%
% [ x(nLags,iCol)   x(nLags-1,iCol) ... x(2,iCol) x(1,iCol) ]
% [ x(nLags+1,iCol) x(nLags,iCol)   ... x(3,iCol) x(2,iCol) ]
% [ ...
% [ x(nPnts,iCol)   x(nPnts-1,iCol) ... x(xe-nLags+2,iCol) x(nPnts-nLags+1,iCol) ]
%

%%
[nPnts,nCol] = size(x);
n = nPnts - nLags + 1;

X = nan(n, nLags * nCol,'like',x);

idx = (1:n);
dimOffset = ((1:nCol)-1)*nLags;

for iLag = 1:nLags
    X(:,dimOffset + iLag) = x(idx - iLag + nLags,:);
end

end
%
%