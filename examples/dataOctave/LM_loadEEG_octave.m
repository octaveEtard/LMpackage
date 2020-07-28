function [data,iB] = LM_loadEEG_octave(c)
% Input c is a cell containing the folder path and file name of an EEG file
% to load.
%
% This function simply loads the file, and return the EEG data into a
% matrix of size [nPnts,nChannels].
%
% It also returns a vector iB containing the indices of stimulus onset for
% all stimuli of interest. This vector should be sorted in the same order
% as the function loading features returns them.
%
% Here this information is already stored
% in the EEG dataset.
%
c = c{1};

EEG = pop_loadset(c{2},c{1});

% index of stimulus onset
eventType = {EEG.event(:).type};
iB = find(strcmp(eventType, 'storyBegin'), 1);

iB = EEG.event(iB).latency;

% data matrix, each column = 1 channel
% make sure channels are returned in the order specified by c{3}
data = LM_reorderMatrix(EEG.data,{EEG.chanlocs(:).labels},c{3},1)';

end
%
%