function saveName = save(data, saveName, saveFolder, overwrite, verbose)
%
% LM.save
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Save 'data' in file 'saveName', in folder 'saveFolder'.
%
%   Usage: proper_save(data, saveName, ...)
%
%   Required input:
%       data        data to be saved: structure or general variable.
%
%   Optional inputs:
%       saveFolder  string  default = '' (current folder)
%       overwrite   logical default = false force overwrite
%       verbose     logical default = false display all operations performed
%
% If 'data' is a structure, and does not include a 'timestamp' field, one
% is added.
%
% If 'data' is not a structure, it will be saved such that the field in the
% resulting structure (.mat) is the name with which proper_save was called.
%
% If 'saveName' doesn't include an extension, it defaults to '.mat'.
%
% If 'saveFolder' does not exist, the necessary folders will be created.
%
% If 'saveName' already exists in 'saveFolder', and 'overwrite' is not set,
% the date (yyyy_mm_dd_HH_MM_SS_) will be prepended to 'saveName', and a
% warning will be displayed.
%
% Setting 'verbose' to true will display a message for the following
% operations: creation of folders, overwriting of files, saving.
%
%
%% Setting default values
% save in current folder
if nargin < 3
    saveFolder = '.';
end
% no overwrite by default
if nargin < 4
    overwrite = 0;
end
% non verbose by default
if nargin < 5
    verbose = 0;
end


%%
[~,saveName,ext] = fileparts(saveName);
% '.mat' if extension not specified
if isempty(ext)
    ext = '.mat';
end
saveName = [saveName, ext];


%% Checking whether saveFolder exists & creating it if needed
% and checking whether saveName already exists, and resolve conflict
% potential warning messages
[saveName,msg] = LM.saveSanityChecks(saveFolder,saveName,overwrite,verbose);


%% Saving
full_path = fullfile(saveFolder, saveName);
if verbose
    fprintf('Saving data in %s\n',full_path);
end
if isstruct(data)
    % timestamping
    if ~isfield(data, 'timestamp')
        data.timestamp = datestr(now,'yyyy_mm_dd_HH_MM_SS');
    end
    save(full_path, '-struct','data','-v7.3');
    fprintf('%s',msg);
else
    % creating a structure so that the variables are saved with their own name
    d.(inputname(1)) = data;
    save(full_path,'-struct', 'd','-v7.3');
    fprintf('%s',msg);
end
end
%
%