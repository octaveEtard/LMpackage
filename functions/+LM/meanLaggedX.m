function mX_row = meanLaggedX(xF,nLags,xb,nPnts,unpad)
%
% meanLaggedX
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Given a timeseries x, represented by its Fourier transform xF, compute
% the mean of the lagged feature matrix X without actually forming it.
%
if unpad
    nFFT = size(xF,1);

    % mean(X,1) can be computed by convolving x with a rectangular window
    % of size n
    mX_row = ifft( xF .* fft(ones(nPnts,1),nFFT,1),'symmetric');
    mX_row = mX_row((1:nLags) + xb + nPnts - 2, : ) / nPnts;
    mX_row = flip(mX_row,1);
    mX_row = mX_row(:)'; % to be equivalent to mean(X,1)
else
    % same mean in each column
    nFeatures = size(xF,2);
    mX_row  = ones(1,nLags*nFeatures);
    
    for iFeature = 1:nFeatures
        mX_row((1:nLags)+nLags*(iFeature-1)) = xF(1,iFeature) / nPnts;
    end  
end
end
%
%