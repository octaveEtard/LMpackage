function [x,nPad_b,nPad_e] = padCCA(x, minLag_x, maxLag_x, minLag, maxLag)
%
% padCCA
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%

% zero pad x
nCol = size(x,2);
[nPad_b, nPad_e] = LM.nPadCCA(minLag_x, maxLag_x, minLag, maxLag);
x = [zeros(nPad_b, nCol) ; x ; zeros(nPad_e, nCol)];

end
%
%