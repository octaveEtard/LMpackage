function [x,nPad_b,nPad_e] = LM_padx(x,minLag,maxLag)
%
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
nFeatures = size(x,2);

[nPad_b,nPad_e] = LM_nPadX(minLag,maxLag);

% zero pad x
x = [zeros(nPad_b,nFeatures) ; x ; zeros(nPad_e,nFeatures)];

end