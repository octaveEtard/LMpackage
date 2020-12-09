function opt = laggedDims(nx,iB,ny,minLag,maxLag)
%
% laggedDims
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Compute relevant indices to form matrices X & y or XtX & Xty.
%
%
% Input:
%
% nx [integer]: number of points in x
% iB [integer]: index of the point in y corresponding to the first point
%                   of x
% ny [integer]: number of points in y
% minLag, maxLag [integers]: extent of the window in x that should be
% considered to model each point in y (minLag <= maxLag is required).
%
%
% # Behaviour:
% - If padding:
%
% Use window in x such that the first window ends at index 1 and the last
% one begins at index nx (full coverage of x / all the points of x are
% used).
%
% If the corresponding points in y do not exist, y is padded with 0 as
% well.
%
% - If no padding is desired:
%
% Use windows such that the first begins at index 1 and the last one ends
% at index nx. This is then corrected to ensure the corresponding points
% are in y.
%
%
% Output (in opt structure):
%
% xb, xe: indices of the first point of the first window / last point of
% the last window that can be used without padding x (i.e.first and last
% points of x that can used with no padding).
%
% yb, ye: indices of the points of y modelled by the first and last window
% with no padding(i.e.first and last point of y that are used).
%
% yb_pad_top, ye_pad_top: indices of the first / last points of y
% corresponding to padded windows at the begining of x
%
% yb_pad_top, ye_pad_top: indices of the first / last point in y
% corresponding to padded windows at the end of x
%
% nZeros_top_y, nZeros_bottom_y: number of zeros to add at the begining
% and end to pad y 
%
%
%% Sanity check
assert(minLag <= maxLag,'minLag <= maxLag is required!');


%% --- Points to use in unpadded X / Y matrices.
%
% In the unpadded matrix X, we attempt to use windows such that the first
% begins at index 1 and the last one ends at index nx (data in x). This is
% then corrected to ensure the corresponding points are in y.
% The first window then corresponds to points index yb (point in y modelled
% by first window / first used point of y) and xb (point in x where the 
% first window begins / first used point of x).
%
% The same is done for the last window (--> xb, yb)
yb = max(1, iB - minLag);
xb = yb - iB + 1 + minLag;

ye = min(ny,iB + nx - 1 - maxLag); 
xe = ye - iB + 1 + maxLag;

%% Points to use in padded X / Y matrices
% 
% With padding, the first window ends at index 1, and the last window
% begins at index nx (full coverage of x).

% (all the points of x are used)
% xb_pad = 1;
% xe_pad = nx;

% first / last point in y corresponding to padded windows at the begining
% of x
yb_pad_top = max(1, iB - maxLag);
ye_pad_top = yb - 1;
% first / last point in y corresponding to padded windows at the end of x
ye_pad_bottom = min(ny, iB + nx -1 - minLag);
yb_pad_bottom = ye + 1;

% number of zeros to add at the begining / end to pad y 
nZeros_top_y = yb_pad_top - (iB - maxLag);
nZeros_bottom_y = (iB + nx -1 - minLag) - ye_pad_bottom;

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

