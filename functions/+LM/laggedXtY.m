function [XtY,x,y] = laggedXtY(x,y,opt,...
    isXFFT,mX,n_mX,Xtop,Xbottom,...
    isYFFT,mY,n_mY,Ytop,Ybottom)
%
% LM.laggedXty
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
%% Preparing input arguments
%
% in this case compute XtX
if isempty(y)
    opt.maxLag_y = opt.maxLag_x;
    opt.minLag_y = opt.minLag_x;
    
    isYFFT = isXFFT;
    auto = true; % auto-correlation
else
    auto = false;
end

% lags
nLags_x = opt.maxLag_x - opt.minLag_x + 1;
nLags_y = opt.maxLag_y - opt.minLag_y + 1;

minLag = min(opt.minLag_x, opt.minLag_y);
maxLag = max(opt.maxLag_x, opt.maxLag_y);
nLags = maxLag - minLag + 1;

% Fourier transform of x and y
if ~(isXFFT || isYFFT)
    nPnts = size(x,1);
    
    if nLags_x == 1 && nLags_y == 1
        error('!'); % FIXME
    elseif nLags_x == 1 || nLags_y == 1
        nPad = max(nLags_x,nLags_y) - 1;
    else
        nPad = nLags - 1;
    end
    nFFT = 2^nextpow2( nPnts + nPad - 1 );
    
elseif isXFFT
    nFFT = size(x,1);
else
    nFFT = size(y,1);
end

if ~isXFFT
    x = fft(x,nFFT,1);
end

if auto
    y = x;
elseif ~isYFFT
    error('!'); % FIXME this is not how fft(y) should be computed
    y = fft(y,nFFT,1);
end

nDims_x = size(x,2);
nDims_y = size(y,2);

dLag = opt.maxLag_y - opt.maxLag_x;

%%
XtY = nan(nDims_x * nLags_x, nDims_y * nLags_y, 'double');

% XtY = X' * Y
% Compute the cross-correlation between x(:,i) and y(:,j) for all (i,j)
% pairs, and fill in XtY

for iDim_x = 1:nDims_x
    % xc x(:,iDim_x) with y(:,1 ... nDims_y)
    xc = ifft(conj(y) .* x(:,iDim_x),nFFT,1,'symmetric');
    % size 2 * nLags - 1 with 0 lag in the middle at index nLags
    xc = [xc( (1:(nLags-1)) - nLags + 1 + nFFT,:); xc(1:nLags,:)];
    % TODO FIXME flip necessary?
    
    if auto
        % making sure XtX is symmetric ; xc(:,1) should be, but there may be
        % numerical differences
        xc(1:nLags,1) = flip(xc(nLags:end,1),1);
    else
        xc = flip(xc,1);
    end
    
    % fill block-rows iFeature
    dimOffset_x = ((1:nDims_x) - 1) * nLags_x;
    dimOffset_y = ((1:nDims_y) - 1) * nLags_y;
    
    for iLag_y = 1:nLags_y
        XtY((1:nLags_x) + dimOffset_x(1),iLag_y + dimOffset_y) = xc( (1:nLags_x) + nLags - 1 - iLag_y + 1 + dLag,:);
    end
end

% TODO merge this with LM.laggedXtX or LM.laggedXty ?

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

if auto
    y = [];
end
end