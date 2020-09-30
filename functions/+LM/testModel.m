function [CC,MSE] = testModel(model,stimOpt,EEGopt,opt,type,mX,mY)
%
% LM.testModel
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
if nargin < 6
    mX = [];
    mY = [];
end

if numel(stimOpt) == 1
    switch type
        case 'forward'
            [CC,MSE] = LM.forward_testModel(model,stimOpt,EEGopt,opt,mX,mY);
        case 'backward'
            [CC,MSE] = LM.backward_testModel(model,stimOpt,EEGopt,opt,mX,mY);
    end
    return;
end

sizeStim = size(stimOpt);
nSub = size(EEGopt,ndims(EEGopt));

if sizeStim(end) == 1
    sizeStim = sizeStim(1:(end-1));
end

assert( all( size(EEGopt) == [sizeStim,nSub] ),...
    'stimOpt & EEGopt dimensions do not match!');

nStimLoad = numel(stimOpt);
nStimPerFile = opt.nStimPerFile;

% number of points to use to measure correlation / MSE
nPntsPerf = opt.nPntsPerf;
nPerfSize = numel(nPntsPerf);

%% Preallocation
CC = cell(nPerfSize,nStimPerFile,nStimLoad);
MSE = cell(nPerfSize,nStimPerFile,nStimLoad);


%%
idxEEGopt = nStimLoad * ((1:nSub)-1);

for iStimulus = 1:nStimLoad
    
    switch type
        case 'forward'
            [CC(:,:,iStimulus),MSE(:,:,iStimulus)] = LM.forward_testModel(model,...
                stimOpt(iStimulus),...
                EEGopt(idxEEGopt + iStimulus),...
                opt,mX,mY);
            
        case 'backward'
            [CC(:,:,iStimulus),MSE(:,:,iStimulus)] = LM.backward_testModel(model,...
                stimOpt(iStimulus),...
                EEGopt(idxEEGopt + iStimulus),...
                opt,mX,mY);
    end
end

CC = reshape(CC,[nPerfSize,nStimPerFile,sizeStim]);
MSE = reshape(MSE,[nPerfSize,nStimPerFile,sizeStim]);

end
%
%