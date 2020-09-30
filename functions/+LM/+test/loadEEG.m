function [data,iB] = loadEEG(c)
%
% LM.test.loadEEG
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Input c is a cell containing the folder path and file name of an EEG file
% to load.
%
% This function simply loads the file, and return the EEG data into a
% matrix of size [nPnts,nChannels].
%
% It also returns a vector iB containing the indices of stimulus onset for
% all stimuli of interest. In this case this information is already stored
% in the EEG dataset.
%
c = c{1};

EEG = load(fullfile(c{1},c{2}));
% data matrix, each column = 1 channel
data = EEG.data;

% indices of stimulus onset
% stimBeginIdx are integers identifying which stimulus was played
iB = EEG.stimBeginIdx;

% returning iB in a specific stimulus order, regardless of presentation
% order

if 2 < numel(c)
    % either:
    % returning specific stimID only in some given order
    [iSort,~] = find( c{3}(:)' == EEG.stimID(:) );
    assert(numel(iSort) == numel(c{3}));
else
    % or:
    % returning all iB, sorted by stimID
    [~,iSort] = sort(EEG.stimID,'ascend');
end
    iB = iB(iSort);
end
%
%