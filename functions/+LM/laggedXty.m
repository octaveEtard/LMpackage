function [Xty,mY,n_mY] = laggedXty(xF,y,minLag,maxLag,...
    mX,Xtop,Xbottom,...
    isYFFT,mY,n_mY,Ytop,Ybottom,...
    opt)
%
% LM.laggedXty
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Compute matrix Xty = X' * y corresponding to a convolutive linear model
% y = X * B with X lagged matrix of the time series x, with lags spanning
% minLag to maxLag (x relative to y).
%
% Xty is the feature - output cross-correlation matrix. Given
% nLags = maxLag - minLag + 1, nFeatures = size(xF,2), and
% nOut = size(y,2), Xty if of size [nLags * nFeatures, nOut], and is
% structured by blocks forming pairwise feature - output cross-correlation
% timeseries.
%
% All the computations are implemented through FFT transforms and without
% explicitely computing X, which is advantageous if the number of
% observations is large.
%
% If opt.removeMean == true, Xty is as if X and y had their column-wise
% means removed before computations (fitting a model without offset).
%
% If opt.unpad.do == true, Xty uses only points within x & y that can be
% included with no padding. Otherwise, all the points of x are used, which
% necessitates zero-padding x with nLags-1 zeros at the edges, and
% potentially zero-padding of y depending on the dimensions.
%
% See also:
%
%   LM.laggedXty to make Xty.
%
%   LM.laggedX & LM.laggedY to make the corresponding X and y matrices.
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
%%
nLags = maxLag - minLag + 1;

iB = opt.iB;

[nFFT,nFeatures] = size(xF);

yb = max(1,iB-maxLag);

if ~isYFFT
    % conjugated of Fourier Transform of y
    [y,mY,n_mY,Ytop,Ybottom] = LM.computeYFFT(y,minLag,maxLag,nFFT,opt);
    %
    % otherwise Ytop, Ybottom, mY and n_mY have to be provided as needed
end

lags = (maxLag:-1:minLag) + (yb - iB);
lags = LM.lagToIndex(lags,nFFT);

nOut = size(y,2);
Xty = nan(nLags*nFeatures,nOut,'double');

for iFeature = 1:nFeatures
    % xc x(:,iFeature) with y
    xc = ifft(xF(:,iFeature) .* y,nFFT,1,'symmetric');
    xc = xc(lags,:);
    
    % fill block-rows iFeature
    dimOffset = (iFeature-1)*nLags;
    Xty((1:nLags) + dimOffset,:) = xc;
end


%%
% Compute Xty as if the data was not padded by removing the top and bottom
% padding

if opt.unpad.do
    Xty = Xty - Xtop' * Ytop - Xbottom' * Ybottom;
end


%% Centering
% These operations are the equivalent of centering X and y before computing
% XtX and Xty:
%
% X = X - mean(X,1);
% y = y - mean(y,1); Xty = X' * y;
%
if opt.removeMean
    Xty = Xty - n_mY * (mX' .* mY);
end

end