
function [eeg,iB] = makeResponse(feature,impResponse,nEdges,nChan,noiseAmp)
%
% LM.test.makeIR
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Convolves feature with impResponse and add some noise. Add nEdges point
% at the begining / end of the generated response ('eeg'). 'iB' is the
% index of the point in 'eeg' corresponding to the first point of 'feature'.
%
% impResponse should be centred on 0.
%
assert(mod(size(impResponse,1),2) == 1);

eeg = [ zeros(nEdges,1) ; ...
    sum(LM.convfft(feature,impResponse,'full'),2) ; ...
    zeros(nEdges,1) ];

nPnts = size(eeg,1);
eeg = eeg + noiseAmp * randn(nPnts,nChan);

iB = (size(impResponse,1)-1) / 2 + 1 + nEdges;

end