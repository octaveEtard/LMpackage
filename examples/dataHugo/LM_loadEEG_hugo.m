function [data,iB] = LM_loadEEG_hugo(c)
% Input c is a cell containing:
%
% c{1}{1} EEG folder path
% c{1}{2} EEG file name
% c{1}{3} order in which to return the channels
% c{1}{4} full path to onset matrix
% c{1}{5} subject index in onset matrix
% c{1}{6} indices of story onset to return
%
% This function loads the EEG file, and return the EEG data into a
% matrix of size [nPnts,nChannels].
%
% It also returns a vector iB containing the indices of stimulus onset for
% all stimuli of interest. This vector should be sorted in the same order
% as the function loading features returns them.
%
% Here this information is obtained from the onset matrix given by c{1}{3} 
%
c = c{1};

% EEG dataset
EEG = pop_loadset(c{2},c{1});

% indices of stimulus onset
d = load(c{4});
iB = ceil( d.onsets(c{5},c{6}) *  EEG.srate) + 1; % from s to samples

% data matrix, each column = 1 channel
% make sure channels are returned in the order specified by c{3}
data = LM_reorderMatrix(EEG.data,{EEG.chanlocs(:).labels},c{3},1)';

end