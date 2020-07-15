function [XtX,Xty] = LM_forward_crossMatrices(stimOpt,EEGopt,opt)
%
% LM_forward_crossMatrices
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
minLag = opt.minLag;
maxLag = opt.maxLag;
nLags = maxLag - minLag + 1;

nChan = opt.nChan;
nStimPerFile = opt.nStimPerFile;
nFeatures = opt.nFeatures;

nSub = numel(EEGopt);

if nStimPerFile == 1
    opt.sumStim = true; % equivalent
end


%% Preallocation
if opt.sumStim
    % XtX will require 8 * (nLags*nFeatures)^2 bytes of memory
    XtX = zeros([nLags*nFeatures,nLags*nFeatures],'double');
    s = [];
else
    % XtX will require 8 * (nLags*nFeatures)^2 * nStimuli  bytes of memory
    XtX = zeros([nLags*nFeatures,nLags*nFeatures,nStimPerFile],'double');
    s = nStimPerFile;
end

if ~opt.sumSub
    s = [s,nSub];
end

Xty = zeros([nLags*nFeatures,nChan,s],'double');


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
mX = zeros(nStimPerFile,nFeatures*nLags);
Xtop = cell(nStimPerFile,1);
Xbottom = cell(nStimPerFile,1);



%% Making Xty for each subject & stimulus
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
            unpad = LM_laggedDims(opt.nx,opt.iB,size(response,1),minLag,maxLag);
            
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
            [XtX_,xF{iStimulus},mX(iStimulus,:),...
                Xtop{iStimulus},Xbottom{iStimulus}] = LM_laggedXtX(feature{iStimulus},minLag,maxLag,opt);
            
            if opt.sumStim
                XtX = XtX + XtX_;
            else
                XtX(:,:,iStimulus) = XtX_;
            end
        end
        
        [Xty_,~] = LM_laggedXty(xF{iStimulus},response,minLag,maxLag,...
            mX(iStimulus,:),Xtop{iStimulus},Xbottom{iStimulus},...
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
end
%
%