function [CC,MSE] = backward_testModel(model,stimOpt,EEGopt,opt,mX,mY)
%
% LM.backward_testModel
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
minLag = opt.minLag;
maxLag = opt.maxLag;
nLags = maxLag - minLag + 1;

nChan = opt.nChan;
nStimPerFile = opt.nStimPerFile;
nFeatures = opt.nFeatures;

% number of points to use to measure correlation / MSE
nPntsPerf = opt.nPntsPerf;

[nOutModel,nReg] = size(model,[2,3]);

[nPad_b,nPad_e] = LM.nPadX(minLag,maxLag);
model = reshape(model,[nLags,nChan,nOutModel,nReg]);

% sanity check
assert( nOutModel == 1 || nOutModel == nFeatures);

nSub = numel(EEGopt);
nPerfSize = numel(nPntsPerf);


%% Preallocation
CC = cell(nPerfSize,nStimPerFile);
MSE = cell(nPerfSize,nStimPerFile);


%% Loading feature representation for all stimuli
% feature sould be a matrix of size [~,nFeatures] or cell with nStimuli
% elements containing matrices of size [~,nFeatures]
feature = opt.getStimulus(stimOpt);

if ~iscell(feature)
    feature = {feature};
end

if opt.removeMean
    mc = nan(1,nChan,nOutModel,nReg);
    mX = reshape(mX,nLags,nChan);
    mc(1,:,:,:) = sum(mX .* model,1);

    mY = reshape(mY,1,nFeatures);
else
    mc = 0;
end

nPnts = nan(nStimPerFile,1);
for iStimulus = 1:nStimPerFile
    % same stimuli for all subjects, computing relevant variables here
    nPnts(iStimulus) = size(feature{iStimulus},1);
    if opt.unpad.do
        n = nPnts(iStimulus);
    else
        feature{iStimulus} = LM.pady(feature{iStimulus},nPnts(iStimulus),nPad_b,nPad_e,1);
        n = nPnts(iStimulus) + 2 * nLags - 2;
    end
    feature{iStimulus} = LM.laggedY(feature{iStimulus},minLag,maxLag,1,n);
    if opt.removeMean
        % feature{iStimulus} = feature{iStimulus}  - mean(feature{iStimulus},1);
        feature{iStimulus} = feature{iStimulus}  - mY;
    end
end


for iSub = 1:nSub
    % response should be a matrix of size [~,nOut]
    % iB should be an array with nStimuli elements containing the index
    % of stimulus onset in response
    [response,iB] = opt.getResponse(EEGopt(iSub));
    
    for iStimulus = 1:nStimPerFile
        nx = nPnts(iStimulus);
        iB_ = iB(iStimulus);
        
        % ------
        % The mean is a learning parameter that should be learnt on the
        % training data if it is to be removed -- otherwise it would
        % require knowledge of all the testing data beforehand, even if
        % testing on smaller slices
        %
        %         if opt.removeMean
        %             mc = nan(1,nChan,nOutModel,nReg);
        %
        %             % compute the mean of corresponding X matrix
        %             if opt.unpad.do
        %                 nFFT = 2^nextpow2( nx + nLags - 1 );
        %                 xF = fft(response((1:nx)+iB_-1,:),nFFT,1);
        %                 n = unpad.xe - unpad.xb - nLags + 2;
        %                 mX = LM.meanLaggedX(xF,nLags,unpad.xb,n,opt.unpad.do);
        %                 mX = reshape(mX,nLags,nChan);
        %             else
        %                 n = nx + nLags - 1;
        %                 mX = repmat(sum(response((1:nx)+iB_-1,:),1) / n,[nLags,1]);
        %             end
        %             mc(1,:,:,:) = sum(mX .* model, 1);
        %         else
        %             mc = 0;
        %         end
        % ------
        
        if opt.unpad.do
            unpad = LM.laggedDims(nx,1,nx,minLag,maxLag);
            pred = squeeze(sum(LM.convfft(response((unpad.xb:unpad.xe)-1+iB_,:),...
                model,'valid') - mc,2));
        else
            pred = squeeze(sum( LM.convfft(response((1:nx)+iB_-1,:),model,'full') - mc,2));
        end
        
        for iPerfSize = 1:nPerfSize
            nPntsWin = nPntsPerf(iPerfSize);
            N = size(pred,1);
            
            if nPntsWin < 0
                nPntsWin = N;
            end
            
            nWin = floor(N / nPntsWin);
            
            if iSub == 1
                % initialisation
                CC{iPerfSize,iStimulus} = nan(nWin,nFeatures,nReg,nSub);
                MSE{iPerfSize,iStimulus} = nan(nWin,nFeatures,nReg,nSub);
            end
            
            for iWin = 1:nWin
                idx = (1:nPntsWin) + (iWin-1) * iWin;
                
                if nOutModel == 1
                    CC{iPerfSize,iStimulus}(iWin,:,:,iSub) = corr(feature{iStimulus}(idx,:),pred(idx,:));
                else
                    for iOut = 1:nOutModel
                        CC{iPerfSize,iStimulus}(iWin,iOut,:,iSub) = corr(feature{iStimulus}(idx,iOut),squeeze(pred(idx,iOut,:)));
                    end
                end
                MSE{iPerfSize,iStimulus}(iWin,:,:,iSub) = shiftdim( sum((feature{iStimulus}(idx,:) - pred(idx,:,:)).^2,1) ./ sum(feature{iStimulus}(idx,:).^2,1),1 );
            end
        end
    end
end
%
%
end
%
%