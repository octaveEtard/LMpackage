function y = pady(y,nx,nPad_b,nPad_e,iB)
%
% pady
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
[ny,nOut] = size(y);
yb = max(1,iB - nPad_b);
ye = min(ny, iB + nx - 1 + nPad_e);

nPad_b = max(0,nPad_b - iB + 1);
nPad_e = max(0,nPad_e - ny + iB + nx - 1);

% pad y
y = [zeros(nPad_b,nOut) ; y(yb:ye,:) ; zeros(nPad_e,nOut)];

end
%
%