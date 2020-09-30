function slices = distribute(n,nSlices)
%
% distribute
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Distribute n in nSlices parts such that sum(slices) == n &
% all element of slices <= n/nSlices
%
k = floor(n/nSlices);

slices = zeros(nSlices,1);
slices(:) = k;
% remainder
r = n - k * nSlices;

if 0 < r
    slices(1:r) = slices(1:r) + 1;
end
%
end