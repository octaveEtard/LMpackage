function features = LM_loadFeature_hugo(filePath)

filePath = filePath{1};
nFiles = numel(filePath);
features = cell(nFiles,1);

for iFile = 1:nFiles
    
    d = load(filePath{iFile});
    features{iFile} = d.story;
    
end

end