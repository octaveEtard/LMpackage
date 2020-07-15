function [x,y] = LM_pad(x,y,minLag,maxLag,iB)
%
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
if nargin < 5
    iB = 1;
end

nx = size(x,1);
[x,nPad_b,nPad_e] = LM_padx(x,minLag,maxLag);

y = LM_pady(y,nx,nPad_b,nPad_e,iB);

end