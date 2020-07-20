function features = LM_loadFeature_hugo(filePath)
% Input filePath is a cell containing paths towards nFiles feature files
% describing each stimulus in one EEG file.
%
% This function simply loads these file, and returns the features.
%
filePath = filePath{1};
nFiles = numel(filePath);
features = cell(nFiles,1);

for iFile = 1:nFiles
    
    d = load(filePath{iFile});
    features{iFile} = d.story;
    
end

end