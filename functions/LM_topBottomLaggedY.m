function [Ytop,Ybottom] =  LM_topBottomLaggedY(y,yb_pad_top,ye_pad_top,...
    yb_pad_bottom,ye_pad_bottom,...
    nPad_top,nPad_bottom)
%
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%

nOut = size(y,2);


nTop = ye_pad_top - yb_pad_top + 1 + nPad_top;
Ytop = zeros(nTop,nOut);
Ytop((nPad_top+1):end,:) = y(yb_pad_top:ye_pad_top,:);


nBottom = ye_pad_bottom - yb_pad_bottom + 1 + nPad_bottom;
Ybottom = zeros(nBottom,nOut);
Ybottom(1:(end-nPad_bottom),:) = y(yb_pad_bottom:ye_pad_bottom,:);

end
%
%