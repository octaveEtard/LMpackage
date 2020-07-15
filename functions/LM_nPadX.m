function [nPad_b,nPad_e] = LM_nPadX(minLag,maxLag)
%
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
nLags = maxLag - minLag + 1;

nPad_b = max(nLags-1, nLags-1+minLag);
nPad_e = max(nLags-1, nLags-1-maxLag);

end