
function [eeg,iB] = LM_testing_makeResponse(feature,impResponse,nEdges,nChan,noiseAmp)

assert(mod(size(impResponse,1),2) == 1);

eeg = [ zeros(nEdges,1) ; ...
    sum(LM_convfft(feature,impResponse,'full'),2) ; ...
    zeros(nEdges,1) ];

nPnts = size(eeg,1);
eeg = eeg + noiseAmp * randn(nPnts,nChan);

iB = (size(impResponse,1)-1) / 2 + 1 + nEdges;

end