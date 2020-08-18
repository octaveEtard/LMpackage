function [XtX,Xty,mX,mY,N] = LM_crossMatrices(stimOpt,EEGopt,opt,type)
%
% LM_crossMatrices
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
%
% if opt.sumSub && 'forward', XtX needs to be multiplied by nSub (same XtX
% accumulated over subjects) for correct normalisation
%
if numel(stimOpt) == 1
    switch type
        case 'forward'
            [XtX,Xty,mX,mY,N] = LM_forward_crossMatrices(stimOpt,EEGopt,opt);
        case 'backward'
            [XtX,Xty,mX,mY,N] = LM_backward_crossMatrices(stimOpt,EEGopt,opt);
    end
    return;
end

minLag = opt.minLag;
maxLag = opt.maxLag;
nLags = maxLag - minLag + 1;

nFeatures = opt.nFeatures;
nChan = opt.nChan;

nStimLoad = numel(stimOpt);
nSub = size(EEGopt,ndims(EEGopt));

sumFrom = opt.sumFrom;


sizeStim = size(stimOpt);
if sizeStim(end) == 1
    sizeStim = sizeStim(1:(end-1));
end

assert( all( size(EEGopt) == [sizeStim,nSub] ),...
    'stimOpt & EEGopt dimensions do not match!');


%% Preallocation
% ---determining size of the ouput
if opt.sumStim || all(opt.nStimPerFile == 1)
    oneStimOutPerFile = true;
    s = [];
else
    % number of stimuli in each stimulus file
    % in this case has to be the same for all elements of stimOpt
    s = opt.nStimPerFile(1);
    oneStimOutPerFile = false;
end

if 0 < sumFrom
    s = [s,sizeStim(1:(sumFrom-1))];
else
    s = [s,sizeStim];
end

if strcmp(type,'forward')
    % --- XtX will require 8 * prod(s) * (nLags*nFeatures)^2 bytes of memory
    XtX = zeros([nLags*nFeatures,nLags*nFeatures,s],'double');
    mX = zeros([nFeatures*nLags,s]);
else
    mY = zeros([nFeatures,s],'double');
end

% number of points (rows) that would be in the equivalent X or Y matrices
% if opt.sumSub, this needs to be multiplied by nSub
N = zeros([1,s]);

if ~opt.sumSub
    s = [s,nSub];
end

switch type
    case 'forward'
        Xty = zeros([nLags*nFeatures,nChan,s],'double');
        mY = zeros([nChan,s]);
        
    case 'backward'
        % --- XtX will require 8 * prod(s) * (nChan*nLags)^2 bytes of memory
        XtX = zeros([nChan*nLags,nChan*nLags,s],'double');
        Xty = zeros([nChan*nLags,nFeatures,s],'double');
        mX = zeros([nChan*nLags,s],'double');
end

if ~oneStimOutPerFile
    s = s(2:end);
end
if numel(s) == 1
    s = [s,1];
end

if 1 < nStimLoad && numel(opt.nStimPerFile) == 1
    % expands
    nStimPerFile = opt.nStimPerFile * ones(nStimLoad,1);
else
    nStimPerFile = opt.nStimPerFile;
end

%%
idxEEGopt = nStimLoad * ((1:nSub)-1);

for iStimulus = 1:nStimLoad
    
    opt.nStimPerFile = nStimPerFile(iStimulus);
    
    % XtX & Xty for all stimuli in stimOpt(iStimulus) and all subjects
    switch type
        case 'forward'
            [XtX_,Xty_,...
                mX_,mY_,...
                N_] = LM_forward_crossMatrices(...
                stimOpt(iStimulus),...
                EEGopt(idxEEGopt + iStimulus),...
                opt);
            
        case 'backward'
            [XtX_,Xty_,...
                mX_,mY_,...
                N_] = LM_backward_crossMatrices(...
                stimOpt(iStimulus),...
                EEGopt(idxEEGopt + iStimulus),...
                opt);
    end
    
    % storing: these do not depend on subject
    [iSlice,idx] = LM_idxSlice(sizeStim,sumFrom,iStimulus);
    
    N(1,iSlice) = N(1,iSlice) + N_;

    if strcmp(type,'forward')
        mX(:,iSlice) = mX(:,iSlice) + N_ .* mX_;
        
        if oneStimOutPerFile
            XtX(:,:,iSlice) = XtX(:,:,iSlice) + XtX_;
        else
            XtX(:,:,:,iSlice) = XtX(:,:,:,iSlice) + XtX_;
        end
    else
        mY(:,iSlice) = mY(:,iSlice) + N_ .* mY_;
    end
    
    % storing: these depend on subject
    if ~opt.sumSub
        if sumFrom == 1
            iSlice = 1:nSub;
        else
            % adding subject index
            tmp = arrayfun(@(i) [idx;i], 1:nSub, 'UniformOutput', false);
            % then converting to linear indices
            iSlice = cellfun(@(c) sub2ind(s,c{:}),tmp);
        end
    end
    
    switch type
        case 'forward'
            mY(:,iSlice) = mY(:,iSlice) + N_ .* mY_;
            
            if oneStimOutPerFile
                Xty(:,:,iSlice) = Xty(:,:,iSlice) + Xty_;
            else
                Xty(:,:,:,iSlice) = Xty(:,:,:,iSlice) + Xty_;
            end
            
        case 'backward'
            mX(:,iSlice) = mX(:,iSlice) +  N_ .* mX_;
            
            if oneStimOutPerFile
                XtX(:,:,iSlice) = XtX(:,:,iSlice) + XtX_;
                Xty(:,:,iSlice) = Xty(:,:,iSlice) + Xty_;
            else
                Xty(:,:,:,iSlice) = Xty(:,:,:,iSlice) + Xty_;
                XtX(:,:,:,iSlice) = XtX(:,:,:,iSlice) + XtX_;
            end
    end
end

mX = mX ./ N;
mY = mY ./ N;
%
%
end
%
%