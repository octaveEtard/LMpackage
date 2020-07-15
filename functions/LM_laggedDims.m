function opt = LM_laggedDims(nx,iB,ny,minLag,maxLag)
%
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Compute relevant indices to form matrices X & y or XtX & Xty
%

% points to use in unpadded X / Y matrices
yb = max(1, iB - minLag);
xb = yb - iB + 1 + minLag;

ye = min(ny,iB + nx - 1 - maxLag); 
xe = ye - iB + 1 + maxLag;

% points to use in padded X / Y matrices
% (all the points of x are used)
% xb_pad = 1;
% xe_pad = nx;
yb_pad_top = max(1, iB - maxLag);
ye_pad_top = yb - 1;

ye_pad_bottom = min(ny, iB + nx -1 - minLag);
yb_pad_bottom = ye + 1;

% number of zeros to add
nZeros_top_y = yb_pad_top - (iB - maxLag);
nZeros_bottom_y = (iB + nx -1 - minLag) - ye_pad_bottom;

% nLags = maxLag - minLag + 1;
% nTop = xb - 2 + nLags;
% nBottom = nLags + nPnts - xe - 1;

opt = struct(...
    'xb',xb,...
    'xe',xe,...
    'yb',yb,...
    'ye',ye,...
    'yb_pad_top',yb_pad_top,...
    'ye_pad_top',ye_pad_top,...
    'yb_pad_bottom',yb_pad_bottom,...
    'ye_pad_bottom',ye_pad_bottom,...
    'nZeros_top_y',nZeros_top_y,...
    'nZeros_bottom_y',nZeros_bottom_y);

end

