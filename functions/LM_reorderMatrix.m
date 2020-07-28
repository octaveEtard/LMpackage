function [M,idx,outputChanOrder] = LM_reorderMatrix(M,inputChanOrder,outputChanOrder,dim)
% reorder matrix M from original order given by inputChanOrder
% to outputChanOrder
%
assert(dim == 1 || dim == 2);

idx = findIndexFullStringInCellArray(inputChanOrder,outputChanOrder);
usedChan = ~cellfun(@isempty,idx);
if any(~usedChan)
    warning('Not all channels found in input!');
end
outputChanOrder = outputChanOrder(usedChan);
idx = [idx{:}];

if isempty(M)
    return;
end
switch dim
    case 1
        M = M(idx,:);
    case 2
        M = M(:,idx);
end
end
%
%