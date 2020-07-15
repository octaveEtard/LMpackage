function [data,iB] = LM_testing_loadEEG(c)
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

EEG = pop_loadset(c{2},c{1});
% data matrix, each column = 1 channel
data = EEG.data';

% indices of stimulus onset
% event code are integers identifying which stimulus was played
eventCode = {EEG.event(:).code};
eventCode = [eventCode{:}];

iB = {EEG.event(:).latency};
iB = [iB{:}];

% returning iB sorted (same stimulus order, regardless of presentation
% order)
[~,iSort] = sort(eventCode,'ascend');
iB = iB(iSort);

end
%
%