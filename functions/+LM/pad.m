function [x,y] = pad(x,y,minLag,maxLag,iB)
%
% pad
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Wrapper to call LM.padx & LM.pady together
%
if nargin < 5
    iB = 1;
end

nx = size(x,1);
[x,nPad_b,nPad_e] = LM.padx(x,minLag,maxLag);

y = LM.pady(y,nx,nPad_b,nPad_e,iB);

end