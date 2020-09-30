function [x,nPad_b,nPad_e] = padx(x,minLag,maxLag)
%
% padx
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
nFeatures = size(x,2);

[nPad_b,nPad_e] = LM.nPadX(minLag,maxLag);

% zero pad x
x = [zeros(nPad_b,nFeatures) ; x ; zeros(nPad_e,nFeatures)];

end