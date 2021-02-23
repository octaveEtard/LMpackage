function [CC,MSE] = forward_testModel(model,stimOpt,EEGopt,opt,mX,mY)
%
% LM.forward_testModel
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
model = reshape(model,[nLags,nFeatures,nOutModel,nReg]);

% sanity check
assert( nOutModel == 1 || nOutModel == nChan);

% how many regularisation coefficients of the model to process at once:
% fewer = less memory needed but slower
% more  = faster but require more memory
switch opt.predBatchSize
    case 'all'
        predBatchSize = nReg;
        nBatch = 1;
    otherwise
        predBatchSize = opt.predBatchSize;
        nBatch = ceil(nReg/ predBatchSize);
end


nPerfSize = numel(nPntsPerf);
nSub = numel(EEGopt);


%% Preallocation
CC = cell(nPerfSize,nStimPerFile);
MSE = cell(nPerfSize,nStimPerFile);


%%
% feature sould be a matrix of size [~,nFeatures] or cell with nStimPerFile
% elements containing matrices of size [~,nFeatures]
feature = opt.getStimulus(stimOpt);

if ~iscell(feature)
    feature = {feature};
end

if opt.removeMean
    mc = nan(1,nFeatures,nOutModel,nReg);
    mX = reshape(mX,nLags,nFeatures);
    if numel(mY) == nChan
        mY = reshape(mY,1,nChan);
    elseif numel(mY) == 1
        % in that case using the same mean for all channels
        mY = mY * ones(1,nChan);
    else
        error('The size of mY does not match');
    end
    mc(1,:,:,:) = sum(mX .* model,1);
end

for iSub = 1:nSub
    % response should be a matrix of size [~,nOut]
    % iB should be an array with nStimuli elements containing the index
    % of stimulus onset in response
    [response,iB] = opt.getResponse(EEGopt(iSub));
    
    for iStimulus = 1:nStimPerFile
        nx = size(feature{iStimulus},1);
        
        % ------
        % The mean is a learning parameter that should be learnt on the
        % training data if it is to be removed -- otherwise it would
        % require knowledge of all the testing data beforehand, even if
        % testing on smaller slices
        %
        %         % FIXME no need to repeat these calculations for all sub
        %         if opt.removeMean
        %             mc = nan(1,nFeatures,nOutModel,nReg);
        %
        %             % compute the mean of corresponding X matrix
        %             if opt.unpad.do
        %                 nFFT = 2^nextpow2( nx + nLags - 1 );
        %                 xF = fft(feature{iStimulus},nFFT,1);
        %                 % assuming opt.unpad.xb == 1 && opt.unpad.xe == nx ;
        %                 % that is, enough response data such that full stimulus
        %                 % lenght can be used given the specified lags
        %                 % n = opt.unpad.xe - opt.unpad.xb - nLags + 2;
        %                 n = nx - 1 - nLags + 2;
        %
        %                 mX = LM.meanLaggedX(xF,nLags,1,n,opt.unpad.do);
        %                 mX = reshape(mX,nLags,nFeatures);
        %             else
        %                 n = nx + nLags - 1;
        %                 mX = repmat(sum(feature{iStimulus},1)/ n,[nLags,1]);
        %             end
        %
        %             mc(1,:,:,:) = sum(mX .* model,1);
        %         end
        % ------
        
        for iBatch = 1:nBatch
            %
            % When a large number of regularisation coefficients (& a large number
            % of output dimensions (e.g. channels)) are used, we may have to work
            % by batch as the predictions for all parameters can't be stored in
            % memory
            
            ibReg = predBatchSize * (iBatch-1) + 1;
            ieReg = min(ibReg + predBatchSize - 1,nReg);
            idxReg = ibReg:ieReg;
            nCurrentBatch = ieReg - ibReg + 1;
            
            % predicted response
            if opt.unpad.do
                pred = LM.convfft(feature{iStimulus},model(:,:,:,idxReg),'valid');
            else
                pred = LM.convfft(feature{iStimulus},model(:,:,:,idxReg),'full');
            end
            if opt.removeMean
                pred = sum(pred - mc(:,:,:,idxReg),2);
            else
                pred = sum(pred,2);
            end
            
            pred = squeeze(pred);
            
            % true response
            if opt.unpad.do
                unpad = LM.laggedDims(nx,iB(iStimulus),size(response,1),minLag,maxLag);
                assert( unpad.xb == 1 && unpad.xe == nx);
                y = response;
                iB_ = iB(iStimulus);
                n = nx;
            else
                y = LM.pady(response,nx,nPad_b,nPad_e,iB(iStimulus));
                iB_ = 1;
                n = size(y,1);
            end
            
            y = LM.laggedY(y,minLag,maxLag,iB_,n);
            
            if opt.removeMean
                % Same as before, mean should be learnt
                % y = y - mean(y,1);
                y = y - mY;
            end
            
            
            for iPerfSize = 1:nPerfSize
                nPntsWin = nPntsPerf(iPerfSize);
                N = size(pred,1);
                
                if nPntsWin < 0
                    nPntsWin = N;
                end
                
                nWin = floor(N / nPntsWin);
                
                if iSub == 1 && iBatch == 1
                    % initialisation
                    CC{iPerfSize,iStimulus} = nan(nWin,nChan,nReg,nSub);
                    MSE{iPerfSize,iStimulus} = nan(nWin,nChan,nReg,nSub);
                end
                
                for iWin = 1:nWin
                    idx = (1:nPntsWin) + (iWin-1) * iWin;
                    
                    if nOutModel == 1
                        CC{iPerfSize,iStimulus}(iWin,:,idxReg,iSub) = corr(y(idx,:),pred(idx,:));
                        for iReg = 1:nCurrentBatch
                            MSE{iPerfSize,iStimulus}(iWin,:,idxReg(iReg),iSub) = shiftdim( sum((y(idx,:) - pred(idx,iReg)).^2,1) ./ sum(y(idx,:).^2,1),1 );
                        end
                    else
                        for iOut = 1:nChan
                            CC{iPerfSize,iStimulus}(iWin,iOut,idxReg,iSub) = corr(y(idx,iOut),squeeze(pred(idx,iOut,:)));
                        end
                        MSE{iPerfSize,iStimulus}(iWin,:,idxReg,iSub) = shiftdim( sum((y(idx,:) - pred(idx,:,:)).^2,1) ./ sum(y(idx,:).^2,1),1 );
                    end
                end
            end
        end
    end
end
%
%
end
%
%