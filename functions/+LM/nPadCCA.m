function [nPad_b,nPad_e] = nPadCCA(minLag_x,maxLag_x,minLag,maxLag)
%
% nPadCCA
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
nPad_b = maxLag - minLag_x;
nPad_e = maxLag_x - minLag;

end