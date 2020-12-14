function [XtY] = laggedXt_2(xF,yF,opt)
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
[nFFT,nDims_x] = size(xF);
[nFFT_,nDims_y] = size(yF);

assert(nFFT == nFFT_);

% TODO here or not?
yF = conj(yF);

nLags_x = opt.maxLag_x - opt.minLag_x + 1;
nLags_y = opt.maxLag_y - opt.minLag_y + 1;

minLag = min(opt.minLag_x, opt.minLag_y);
maxLag = max(opt.maxLag_x, opt.maxLag_y);
nLags = maxLag - minLag + 1;

dLag = opt.maxLag_y - opt.maxLag_x;

XtY = nan(nDims_x * nLags_x, nDims_y * nLags_y, 'double');

% XtY = X' * Y
% Compute the cross-correlation between x(:,i) and y(:,j) for all (i,j)
% pairs, and fill in XtY

for iDim_x = 1:nDims_x
    % xc x(:,iDim_x) with y(:,1 ... nDims_y)
    xc = ifft(yF .* xF(:,iDim_x),nFFT,1,'symmetric');
    % size 2 * nLags - 1 with 0 lag in the middle at index nLags
    xc = flip([xc( (1:(nLags-1)) - nLags + 1 + nFFT,:); xc(1:nLags,:)],1);
    
    % fill block-rows iFeature
    dimOffset_x = ((1:nDims_x) - 1) * nLags_x;
    dimOffset_y = ((1:nDims_y) - 1) * nLags_y;
    
    for iLag_y = 1:nLags_y
        XtY((1:nLags_x) + dimOffset_x(1),iLag_y + dimOffset_y) = xc( (1:nLags_x) + nLags - 1 - iLag_y + 1 + dLag,:);
    end
end

% TODO merge this with LM.laggedXtX ?

%%
% Compute Xty as if the data was not padded by removing the top and bottom
% padding
if opt.unpad.do
    XtY = XtY - Xtop' * Ytop - Xbottom' * Ybottom;
end


%% Centering
% These operations are the equivalent of centering X and y before computing
% XtX and Xty:
%
% X = X - mean(X,1);
% y = y - mean(y,1); Xty = X' * y;
%
if opt.removeMean
    XtY = XtY - n_mY * (mX' .* mY);
end

end