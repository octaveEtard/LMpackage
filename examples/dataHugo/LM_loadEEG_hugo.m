function [EEG,iB] = LM_loadEEG_hugo(c)
% Input c is a cell containing:
%
% c{1}{1} EEG folder path
% c{1}{2} EEG file name
% c{1}{3} full path to onset matrix
% c{1}{4} subject index in onset matrix
% c{1}{5} indices of story onset to return
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
%
EEG = pop_loadset(c{2},c{1});

% removing 'Sound' channel
iSound = find(strcmp({EEG.chanlocs(:).labels},'Sound'),1);
iEEG = 1:EEG.nbchan;
iEEG(iSound) = []; % actual EEG channels;


d = load(c{3});
iB = ceil( d.onsets(c{4},c{5}) *  EEG.srate) + 1; % from s to samples

% data matrix, each column = 1 channel
EEG = EEG.data(iEEG,:)';

end