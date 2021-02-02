function lags = lagToIndex(lags,nFFT)
%
% LM.lagToIndex
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Convert lags (e.g. -2,-1,0,1) to the index they are at in the result of
% cross-correlation by FFT.
%
% ! this function does not check that the requested lags are valid
% (abs() < nFFT /2 +1)

% lag 0 is at index 1 followed by positive lags
pLags = 0 <= lags;
lags( pLags ) = lags( pLags ) + 1;

% negative lags are at the end of the array
lags( ~pLags ) = lags( ~pLags ) + nFFT + 1;

end
%
%