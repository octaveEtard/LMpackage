function [XtX,Xty,mX,mY,N,yty] = CCA_crossMatrices(stimOpt,EEGopt,opt)
%
% LM.CCA_crossMatrices
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% CCA matrices with X: stimulus ; y: response
%
minLag_s = opt.stimulus.minLag;
maxLag_s = opt.stimulus.maxLag;
nLags_s = maxLag_s - minLag_s + 1;

minLag_r = opt.response.minLag;
maxLag_r = opt.response.maxLag;
nLags_r = maxLag_r - minLag_r + 1;

nChan = opt.nChan;
nStimPerFile = opt.nStimPerFile;
nFeatures = opt.nFeatures;

nSub = numel(EEGopt);

if nStimPerFile == 1
    opt.sumStim = true; % equivalent
end

% ?
% unpad = opt.unpad.do;
% ?


%% Preallocation
if opt.sumStim
    % XtX will require 8 * (nLags_s * nFeatures)^2 bytes of memory
    XtX = zeros(nLags_s*nFeatures*[1,1],'double');
    s = [];
else
    % XtX will require 8 * (nLags_s*nFeatures)^2 * nStimuli  bytes of memory
    XtX = zeros([nLags_s*nFeatures,nLags_s*nFeatures,nStimPerFile],'double');
    s = nStimPerFile;
end
if ~opt.sumSub
    s = [s,nSub];
end

% yty will require 8 * (nLags_r * nChan)^2 (* nStimPerFile) (* nSub) bytes of memory
yty = zeros([nLags_r * nChan, nLags_r * nChan,s ],'double');

% Xty will require 8 * (nLags_s * nFeatures) * (nLags_r * nChan) (* nStimPerFile) (* nSub) bytes of memory
Xty = zeros([nLags_s * nFeatures, nLags_r * nChan, s],'double');


%% Loading feature representation for all stimuli
% feature sould be a matrix of size [~,nFeatures] or cell with nStimPerFile
% elements containing matrices of size [~,nFeatures]
feature = opt.getStimulus(stimOpt);

if ~iscell(feature)
    feature = {feature};
end


%% Making XtX for each stimulus (for a given stimulus XtX is the same for
% all subjects)
xF = cell(nStimPerFile,1);

% number of points (rows) that would be in the equivalent X or Y matrices
% if opt.sumSub, this needs to be multiplied by nSub
N = nan(1,nStimPerFile);

mX = nan(nFeatures*nLags_s,nStimPerFile,'double');
mY = nan(nLags_r*nChan,nStimPerFile,nSub,'double');

Xtop = cell(nStimPerFile,1);
Xbottom = cell(nStimPerFile,1);


%% Making Xty & yty for each subject & stimulus
xbe = nan(2,nStimPerFile);

for iSub = 1:nSub
    % response should be a matrix of size [~,nOut]
    % iB should be an array with nStimuli elements containing the index
    % of stimulus onset in response
    [response,iB] = opt.getResponse(EEGopt(iSub));
    
    for iStimulus = 1:nStimPerFile
        
        opt.nx = size(feature{iStimulus},1);
        opt.iB = iB(iStimulus);
        
        
        if opt.unpad.do
            unpad = LM.laggedDims(opt.nx,opt.iB,size(response,1),minLag,maxLag);
            
            if iSub == 1
                xbe(1,iStimulus) = unpad.xb;
                xbe(2,iStimulus) = unpad.xe;
            else
                % otherwise the Xtop & Xbottom needed to unpad won't be the
                % same across subjects
                assert( xbe(1,iStimulus) == unpad.xb && ...
                    xbe(2,iStimulus) == unpad.xe );
            end
            
            opt.unpad = unpad;
            opt.unpad.do = true;
        end
        
        if iSub == 1
            % same for all subject, but need unpad structure to be
            % initialised if unpadding
            [XtX_,xF{iStimulus},mX(:,iStimulus),...
                Xtop{iStimulus},Xbottom{iStimulus},N(iStimulus)] = LM.laggedXtX(feature{iStimulus},minLag,maxLag,opt);
            
            if opt.sumStim
                XtX = XtX + XtX_;
            else
                XtX(:,:,iStimulus) = XtX_;
            end
        end
        
        [Xty_,mY(:,iStimulus,iSub)] = LM.laggedXty(xF{iStimulus},response,minLag,maxLag,...
            mX(:,iStimulus)',Xtop{iStimulus},Xbottom{iStimulus},...
            false,[],[],[],[],opt);
        
        if opt.sumSub
            if opt.sumStim
                Xty = Xty + Xty_;
            else
                Xty(:,:,iStimulus) = Xty(:,:,iStimulus) + Xty_;
            end
        else
            if opt.sumStim
                Xty(:,:,iSub) = Xty(:,:,iSub)  + Xty_;
            else
                Xty(:,:,iStimulus,iSub) = Xty(:,:,iStimulus,iSub) + Xty_;
            end
        end
    end
end
%
%
if opt.sumStim
    mX = sum(N .* mX,2) ./ sum(N,2);
    
    if opt.sumSub
        mY = sum(N .* mY,[2,3]) ./ (sum(N,2)*nSub);
    else
        mY = sum(N .* mY,2) ./ sum(N,2);
        mY = permute(mY,[1,3,2]); % drop 2nd dim only
    end
    
    N = sum(N,2);
    
else
    if opt.sumSub
        mY = sum(mY,3) / nSub;
    end
end
%
%
end
%
%