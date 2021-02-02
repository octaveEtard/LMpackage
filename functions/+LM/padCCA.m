function [x,y] = padCCA(x, y, iB, minLag_x, maxLag_x, minLag_y, maxLag_y)
%
% padCCA
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
opt = LM.CCAdims(size(x,1), size(y,1), iB, minLag_x, maxLag_x, minLag_y, maxLag_y);

x = zeroPadMatrix(x,opt.x.nPad_b,opt.x.nPad_e,opt.xb,opt.xe);
y = zeroPadMatrix(y,opt.y.nPad_b,opt.y.nPad_e,opt.yb,opt.ye);
end
%
function x = zeroPadMatrix(x,nPad_b,nPad_e,xb,xe)

nCol = size(x,2);
x = [zeros(nPad_b, nCol) ; x(xb:xe,:) ; zeros(nPad_e, nCol)];

end
%
%