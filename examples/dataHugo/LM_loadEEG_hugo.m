function [EEG,iB] = LM_loadEEG_hugo(c)
%
c = c{1};
% c{1} EEG folder path
% c{2} EEG file name
% c{3} full path to onset matrix
% c{4} subject index in onset matrix
% c{5} indices of story onset to return
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