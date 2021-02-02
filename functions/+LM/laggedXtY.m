function XtY = laggedXtY(xF,yF,minLag_x,maxLag_x,minLag_y,maxLag_y,n_m,...
    mX,Xtop,Xbottom,...
    mY,Ytop,Ybottom,...
    opt)
%
% LM.laggedXty ---- Work in progress -- incomplete code!
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
%% Preparing input arguments
%
% lags
nLags_x = maxLag_x - minLag_x + 1;
nLags_y = maxLag_y - minLag_y + 1;

% minLag = min(minLag_x, minLag_y);
% maxLag = max(maxLag_x, maxLag_y);
% nLags = maxLag - minLag + 1;

[nFFT,nDims_x] = size(xF);
nDims_y = size(yF,2); % size(yF,1) sould be == nFFT too

yF = conj(yF);


%%
lags = ((maxLag_x - minLag_y):-1:(minLag_x - maxLag_y)) + (opt.dims.yb - opt.dims.y.iB);
lags = LM.lagToIndex(lags,nFFT);

% XtY = X' * Y
% Compute the cross-correlation between x(:,i) and y(:,j) for all (i,j)
% pairs, and fill in XtY
dimOffset_y = ((1:nDims_y) - 1) * nLags_y;
XtY = nan(nDims_x * nLags_x, nDims_y * nLags_y, 'double');

for iDim_x = 1:nDims_x
    % xc x(:,iDim_x) with y(:,1 ... nDims_y)
    % y already conj
    xc = ifft(yF .* xF(:,iDim_x),nFFT,1,'symmetric');
    xc = xc(lags,:);
    
    % fill block-rows iFeature
    dimOffset_x = (iDim_x - 1) * nLags_x;
    
    for iLag_y = 1:nLags_y
        % XtY((1:nLags_x) + dimOffset_x,iLag_y + dimOffset_y) = xc( (1:nLags_x) + nLags - 1 - iLag_y + 1 + dLag,:);
        XtY((1:nLags_x) + dimOffset_x,iLag_y + dimOffset_y) = xc( (1:nLags_x) + nLags_y - iLag_y,:);
    end
end


%%
% Compute Xty as if the data was not padded by removing the top and bottom
% padding
if opt.unpad.do
    XtY = XtY - Xtop' * Ytop - Xbottom' * Ybottom;
end


%% Centering
% These operations are the equivalent of centering the lagged matrices X &
% Y before computing XtY:
%
% X = X - mean(X,1);
% Y = y - mean(Y,1);
% XtY = X' * Y;
%
if opt.removeMean
    XtY = XtY - n_m * (mX' .* mY);
end
end
%
%