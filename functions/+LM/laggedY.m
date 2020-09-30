function y = laggedY(y,minLag,maxLag,iB,nx)
%
% laggedY
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Given a timeseries y with [nPnts,nOut] = size(y),
% minLag <= maxLag integers, nx the number of points in a timeseries x, and
% iB the time index to synchronize x and y (i.e. such that y(iB,:)
% corresponds to x(1,:)), make matrix y_ corresponding to a convolutive
% linear model y_ = X * B where each point y_(i,iOut) is obtained as the
% weighted sum over points in x in  the range
% x(i+minLag,iFeature):x(i+maxLag,iFeature) and over each input
% feature (column) in x.
%
% This function does NOT pad the data, and only valid points in x & y are
% considered ; i.e. the matrix y_ is formed using only points in the
% timeseries y, and such that the corresponding matrix X can be formed
% using only points within x.
%
% If padding is desired, it can be achieved by appropriately padding the
% inputs x & y. See 'LM_pad' for example.
%
% With yb and ye the indices of the first and last usable points in y,
% nLags = maxLag - minLag + 1, n = ye - yb  + 1, return lagged
% matrix y_ of size [n,nOut]:
%
% y_ = y(yb:ye,:)
%
% See also 'laggedX which computes the corresponding matrix X, and
% 'laggedXty' which efficiently computes Xty = X' * y_ through FFT
% transforms and without explicitely computing X or y_, which is
% advantageous if the number of observations is large.
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
% ------
%%
ny = size(y,1);

if nargin < 4
    iB = 1;
end
if nargin < 5
    nx = ny;
end

opt = LM.laggedDims(nx,iB,ny,minLag,maxLag);
y = y(opt.yb:opt.ye,:);

end
%
%