function [XtX,xF,mX,Xtop,Xbottom,n_mX] = laggedXtX(x,minLag,maxLag,opt)
%
% LM.laggedXtX
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Compute matrix XtX = X' * X corresponding to a convolutive linear model
% y = X * B with X lagged matrix of the time series x, with lags spanning
% minLag to maxLag (x relative to y).
%
% XtX is the feature auto-correlation matrix. Given
% nLags = maxLag - minLag + 1, and nFeatures = size(x,2), XtX if of size
% [nLags * nFeatures, nLags * nFeatures], and is structured by blocks
% forming pairwise feature cross-correlation matrices.
%
% All the computations are implemented through FFT transforms and without
% explicitely computing X, which is advantageous if the number of
% observations is large.
%
% If opt.removeMean == true, XtX is as if X had its column-wise mean
% removed before computations (fitting a model without offset).
%
% If opt.unpad.do == true, X uses only points within x that can be included
% with no padding. Otherwise, all the points of x are used, which
% necessitates zero-padding x with nLags-1 zeros at the edges.
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
[nPnts,nFeatures] = size(x);
nLags = maxLag - minLag + 1;
% explicit padding value
if isfield(opt,'nPad')
    nPad = opt.nPad;
else
    nPad = nLags - 1; % default value (non-periodic, zero padded signal)
end

XtX = nan(nLags*nFeatures,nLags*nFeatures,'like',x);

nFFT = 2^nextpow2( nPnts + nPad );
xF = fft(x,nFFT,1);

% XtX = X' * X
% Compute the cross-correlation between x(:,i) and x(:,j) for all i <= j
% pairs, and fill in XtX

for iFeature = 1:nFeatures
    % xc x(:,iFeature) with x(:,iFeature ... nFeature)
    xc = ifft(conj(xF(:,iFeature)) .* xF(:,iFeature:end),nFFT,1,'symmetric');
    % size 2 * nLags - 1 with 0 lag in the middle at index nLags
    xc = [xc( (1:(nLags-1)) - nLags + 1 + nFFT,:); xc(1:nLags,:)];
    % making sure XtX is symmetric ; xc(:,1) should be, but there may be
    % numerical differences
    xc(1:nLags,1) = flip(xc(nLags:end,1),1);
    
    % fill block-rows iFeature
    dimOffset = ((iFeature:nFeatures)-1)*nLags;
    for iLag = 1:nLags
        XtX((1:nLags) + dimOffset(1),iLag + dimOffset) = xc( (1:nLags) - 1 + nLags  - iLag + 1,:);
    end
    
    if iFeature < nFeatures
        % fill the symmetric block-colums iFeature
        XtX((iFeature*nLags+1):end, (1:nLags) + dimOffset(1)) = XtX((1:nLags) + dimOffset(1),(iFeature*nLags+1):end)';
    end
end


%%
n_mX = nPnts + nPad;

if opt.unpad.do
    % Compute XtX as if the data was not padded by removing the top and bottom
    % padding
    xb = opt.unpad.xb;
    [Xtop,Xbottom] =  LM.topBottomLaggedX(x,xb,opt.unpad.xe,minLag,maxLag);
    XtX = XtX - (Xtop' * Xtop) - (Xbottom' * Xbottom);
    n_mX = n_mX - size(Xtop,1) - size(Xbottom,1);
else
    Xtop = []; Xbottom = [];
    xb = [];
end


%% Centering
% These operations are the equivalent of centering X before computing XtX
%
% X = X - mean(X,1); XtX = X' * X;
%
if opt.removeMean || 2 < nargout
    mX = LM.meanLaggedX(xF,nLags,xb,n_mX,opt.unpad.do);
end

if opt.removeMean
    XtX = XtX - n_mX * (mX' * mX);
end

end
%
%