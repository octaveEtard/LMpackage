function [nPad_b,nPad_e] = nPadX(minLag,maxLag)
%
% nPadX
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%

% nLags = maxLag - minLag + 1;

% nPad_b = max(nLags-1, nLags-1+minLag);
nPad_b = max(maxLag - minLag, maxLag );  % same

% nPad_e = max(nLags-1, nLags-1-maxLag);
nPad_e = max(maxLag - minLag, -minLag );  % same
end
%
%