function feature = LM_loadFeature_octave(c)
% Input c is a cell containing the path towards one feature file
% describing the stimulus
%
% This function simply loads the file, and returns the feature.
%
d = load(c{1});
feature = d.attended;

end
%
%