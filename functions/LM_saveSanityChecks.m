function [saveName,msg] = LM_saveSanityChecks(saveFolder,saveName,overwrite,verbose)
%
% no overwrite by default
if nargin < 3
    overwrite = 0;
end
% non verbose by default
if nargin < 4
    verbose = 0;
end


%% Checking whether saveFolder exists
full_path = fullfile(saveFolder, saveName);
if ~strcmp(saveFolder,'')
    if ~exist(saveFolder, 'dir')
        if verbose
            fprintf(['Folder ', saveFolder, ' does not exist, creating it.\n']);
        end
        mkdir(saveFolder);
    end
end


%% Checking whether saveName already exists
% potential warning messages
msg = '';
if exist(full_path, 'file')
    if isempty(saveFolder)
        saveFolderMsg = pwd();
    else
        saveFolderMsg = saveFolder;
    end
    if overwrite
        if verbose
            msg = sprintf('%s\n', ['File ', saveName, ' already exists in ', saveFolderMsg, ' - OVERWRITE!']);
        end
    else
        new_name = [datestr(now,'yyyy_mm_dd_HH_MM_SS_'), saveName];
        msg = warning('%s\n', ['File ', saveName, ' already exists in ', saveFolderMsg, ' - New name: ', new_name]);
        saveName = new_name;
    end
end

end