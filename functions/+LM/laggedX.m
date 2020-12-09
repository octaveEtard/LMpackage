function X = laggedX(x,minLag,maxLag,iB,ny)
%
% laggedX
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Given a timeseries x with [nPnts,nFeatures] = size(x),
% minLag <= maxLag integers, ny the number of points in a timeseries y, and
% iB the time index to synchronize x and y (i.e. such that y(iB,:)
% corresponds to x(1,:)), make matrix X corresponding to a convolutive
% linear model y = X * B where each point y(i,iOut) is obtained as the
% weighted sum over points in x in  the range
% x(i+minLag,iFeature):x(i+maxLag,iFeature) and over each input
% feature (column) in x.
%
% This function does NOT pad the data, and only valid points in x & y are
% considered ; i.e. the matrix X is formed using only points in the
% timeseries x, and such that the corresponding matrix y can be formed
% using only points within y.
%
% If padding is desired, it can be achieved by appropriately padding the
% inputs x & y. See 'LM_pad' for example.
%
% With xb and xe the indices of the first and last usable points in x,
% nLags = maxLag - minLag + 1, n = xe - xb - nLags + 2, return lagged
% matrix  X of size [n,nLags * nFeatures]:
%
% X(:,(1:nLags) + (iFeature-1) * nLags) =
%
% [ x(xb+nLags-1,iDim) x(xb+nLags-2,iDim) ... x(xb+1,iDim)   x(xb,iDim)   ]
% [ x(xb+nLags,iDim)   x(xb+nLags-1,iDim) ... x(xb,iDim)     x(xb+1,iDim) ]
% [ ...
% [ x(xe,iDim)   x(xe-1,iDim)  ...  x(xe-nLags+2,iDim) x(xe-nLags+1,iDim) ]
%
% See also 'LM.laggedY which computes the corresponding matrix Y, and
% 'LM.laggedXtX' which efficiently computes XtX = X' * X through FFT
% transforms and without explicitely computing X, which is advantageous if
% the number of observations is large.
%
% ------ Example: -------
% N = 100;
% % x and y some time series
% x = randn(N,3);
% y = randn(N,4);
% 
% iB = 1;
% minLag = -3;
% maxLag = 2;
% 
% % --- Computing XtX and Xty directly:
% opt = struct();
% opt.unpad.do = false;
% opt.removeMean = true;
% 
% [XtX_FFT,xF,mX] = LM.laggedXtX(x,minLag,maxLag,opt);
% 
% opt.iB = iB;
% opt.nx = size(x,1);
% 
% Xty_FFT = LM.laggedXty(xF,y,minLag,maxLag,mX,[],[],false,[],[],[],[],opt);
% 
% % --- Same as:
% % Padding:
% [x,y] = LM.pad(x,y,minLag,maxLag);
% 
% X = LM.laggedX(x,minLag,maxLag,iB,size(y,1));
% Y = LM.laggedY(y,minLag,maxLag,iB,size(x,1));
% 
% X = X - mean(X,1);
% y = y - mean(y,1);
% 
% XtX = X' * X;
% Xty = X' * Y;
% 
% % --- Check:
% max(abs(XtX-XtX_FFT),[],'all')
% max(abs(Xty-Xty_FFT),[],'all')
%
%------
%%
nPnts = size(x,1);

if nargin < 4
    iB = 1;
end
if nargin < 5
    ny = nPnts;
end

opt = LM.laggedDims(nPnts,iB,ny,minLag,maxLag);
% xb = opt.xb;
% xe = opt.xe;

X = LM.lagMatrix(x(opt.xb:opt.xe,:), nLags);

% nLags = maxLag - minLag + 1;
% n = xe - xb - nLags + 2;
% 
% X = nan(n, nLags * nFeatures,'like',x);
% 
% idx = (1:n) + xb - 1;
% dimOffset = ((1:nFeatures)-1)*nLags;
% 
% for iLag = 1:nLags
%     X(:,dimOffset + iLag) = x(idx - iLag + nLags,:);
% end

end
%
%