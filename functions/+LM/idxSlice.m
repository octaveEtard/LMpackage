function [iSlice,idx] = idxSlice(sizeVec,sumFrom,iLinear)
%
% LM.idxSlice
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Given a linear index iLinear in a matrix of size sizeVec, returns the
% linear index iSlice and indices idx corresponding by dropping dimensions
% from index sumFrom on.

if sumFrom == 1
    % in this case all dimensions are dropped
    iSlice = 1;
    idx = {1};
    return
end

% indices corresponding to iLinear
% (we need the cell trick to deal with the varargout & varargin situation)
idx = cell(numel(sizeVec),1);
[idx{:}] = ind2sub(sizeVec,iLinear);

if sumFrom <= 0
    % nothing happens in this case
    %     idx = [idx{:}];
    iSlice = iLinear;
    return;
end

% remaning indices & dimensions
idx = idx(1:(sumFrom-1));
sizeVec = sizeVec(1:(sumFrom-1));

% Matlab wants at least 2 dimensions
if numel(sizeVec) == 1
    sizeVec = [sizeVec,1];
end

iSlice = sub2ind(sizeVec,idx{:});
% idx = [idx{:}];

end
%
%