function indices = findIndexFullStringInCellArray(cellArray, targetStrings)
% findIndexFullStringInCellArray
%
% Input:
% 'cellArray': cell array of strings to search
% 'targetStrings': cell array of strings to look for
%
% Ouput:
% 'indices': indices of cells in 'cellArray' where 'targetString' was found
% cell array of integers, with empty cells where the target string was not
% found
%
% Example: findIndexFullStringInCellArray({'aba', 'bcb', 'a'}, {'a'})
% returns:
% 'indices': {3}
%

function indices = myFind(targetString)
indices = find(strcmp(cellArray,targetString));
end

% myFind = @(targetString) find(strcmp(cellArray,targetString));

% UniformOutput in case one string is not found, or found several times
indices = cellfun(@myFind,targetStrings,'UniformOutput',false);

end

