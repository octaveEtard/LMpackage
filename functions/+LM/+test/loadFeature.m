function feature = loadFeature(c)
%
% LM.test.loadFeature
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Input c is a cell containing n paths towards feature files describing the
% stimuli.
%
% This function simply loads each file, and put the feature data (matrix)
% in a new cell array.
%
c = c{1};
n = numel(c);
feature = cell(n,1);

for iFile = 1:n
    d = load(c{iFile});
    feature{iFile} = d.feature;  
end
end
%
%