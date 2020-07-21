function chanLocs = LM_example_loadChanLocs(nbChan)

if nargin < 1
    nbChan = 64;
end

% assuming this is run from LMpackage/examples/someExample
chanLocs = load(fullfile('..',sprintf('chanLocs-%i.mat',nbChan)));
chanLocs = chanLocs.chanLocs;

end