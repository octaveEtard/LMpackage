% Testing function part of the Linear Model (LM) package.
% Author: Octave Etard
%
% run 'generate_testing_data' to generate synthetic data first
% edit folder name
baseDataFolder = 'testing_data_someTimestamp';

maxDev = zeros(4,1);

nStimPerFile = 4;

minLagT = -500e-3;
maxLagT = 500e-3;

nChan = 5;

Fs = 100;
durStim = 3 * 60;

iCond = 1;
iPart = 1;
iSub = 1;

stimOpt = cell(nStimPerFile,1);
for iStim = 1:nStimPerFile
    stimOpt{iStim} = fullfile(baseDataFolder,'features',...
        sprintf('feature_cond_%i_part_%i_stim_%i.mat',iCond,iPart,iStim));
end
stimOpt = {stimOpt};

EEGopt = {{fullfile(baseDataFolder,'EEG'),...
    sprintf('sub_%i_cond_%i_part_%i.set',iSub,iCond,iPart)}};


opt = struct();
opt.nStimPerFile = nStimPerFile;
opt.getStimulus = @LM_testing_loadFeature;
opt.getResponse = @LM_testing_loadEEG;

opt.sumSub = false;
opt.sumStim = false;

opt.nChan = nChan;
opt.nFeatures = 1;
opt.minLag = floor(minLagT * Fs);
opt.maxLag = ceil(maxLagT * Fs);

opt.nPntsPerf = -1;

trainOpt = struct();
trainOpt.printOut = false;
trainOpt.accumulate = true;
trainOpt.method.name = 'ridge-eig-XtX';
trainOpt.method.lambda = 10.^(-6:0.5:6);
trainOpt.method.normaliseLambda = true;
% alternatively the following 2 lines could be skipped to use default tol
trainOpt.method.removeEig.type = 'tol';
trainOpt.method.removeEig.tol = 0; % default tolerance

for padded = [true,false]
    opt.unpad.do = ~padded;
    
    for removeMean = [true,false]
        
        opt.removeMean = removeMean;
        [XtX,Xty] = LM_crossMatrices(stimOpt,EEGopt,opt,'backward');
        % ------
        
        % ------
        feature = LM_testing_loadFeature(stimOpt);
        [EEG,iB] = LM_testing_loadEEG(EEGopt);
        
        for iStimulus = 1:opt.nStimPerFile
            n = size(feature{iStimulus},1);
            iB_ = iB(iStimulus);
            if padded
                [xp,yp] = LM_pad(EEG((1:n)-1+iB_,:),feature{iStimulus},opt.minLag,opt.maxLag,1);
            else
                xp = EEG((1:n)-1+iB_,:);
                yp = feature{iStimulus};
            end
            
            X = LM_laggedX(xp,opt.minLag,opt.maxLag,1,size(yp,1));
            Y = LM_laggedY(yp,opt.minLag,opt.maxLag,1,size(xp,1));
            
            if removeMean
                X = X - mean(X,1);
                Y = Y - mean(Y,1);
            end
            % ------
            
            % --- these should be the same
            XtX_ = X'*X;
            Xty_ = X'*Y;
            
            maxDev(1) = max(maxDev(1),max(abs((XtX(:,:,iStimulus,iSub)-XtX_) ./ XtX_),[],'all'));
            maxDev(2) = max(maxDev(2),max(abs((Xty(:,:,iStimulus,iSub)-Xty_) ./ Xty_),[],'all'));
            % ---
        end
        % ------
        XtX = sum(XtX,3);
        Xty = sum(Xty,3);
        
        model = LM_fitLinearModel(XtX,Xty,trainOpt);
        [nOutModel,nLambda] = size(model.coeffs,[2,3]);
        
        %         for predBatchSize = [1,ceil(nLambda/3),nLambda]
        predBatchSize = nLambda;
        opt.predBatchSize = predBatchSize;
        [CC,MSE] = LM_testModel(model.coeffs,stimOpt,EEGopt,opt,'backward');
        
        CC = squeeze(permute(vertcat(CC{:}),[2,3,1]));
        MSE = squeeze(permute(vertcat(MSE{:}),[2,3,1]));
        
        % --- prediction with the X and y matrices computed ---
        CC_ = nan(opt.nFeatures,nLambda,opt.nStimPerFile);
        MSE_ = nan(opt.nFeatures,nLambda,opt.nStimPerFile);
        
        for iStimulus = 1:opt.nStimPerFile
            n = size(feature{iStimulus},1);
            iB_ = iB(iStimulus);
            if padded
                [xp,yp] = LM_pad(EEG((1:n)-1+iB_,:),feature{iStimulus},opt.minLag,opt.maxLag,1);
            else
                xp = EEG((1:n)-1+iB_,:);
                yp = feature{iStimulus};
            end
            
            X = LM_laggedX(xp,opt.minLag,opt.maxLag,1,size(yp,1));
            Y = LM_laggedY(yp,opt.minLag,opt.maxLag,1,size(xp,1));
            
            if removeMean
                X = X - mean(X,1);
                Y = Y - mean(Y,1);
            end
            
            for iLambda = 1:nLambda
                pred = X * model.coeffs(:,:,iLambda);
                
                if nOutModel == 1
                    CC_(1,iLambda,iStimulus) = corr(Y,pred);
                    MSE_(1,iLambda,iStimulus) = sum( (Y-pred).^2, 1) ./ sum(Y.^2,1);
                else
                    for iOut = 1:nOutModel
                        CC_(iOut,iLambda,iStimulus) = corr(Y(:,iOut),pred(:,iOut));
                        MSE_(iOut,iLambda,iStimulus) = sum( (Y(:,iOut)-pred(:,iOut)).^2, 1) / sum(Y(:,iOut).^2,1);
                    end
                end
            end
        end
        CC_ = squeeze(CC_);
        MSE_ = squeeze(MSE_);
        % ---
        maxDev(3) = max(maxDev(3),max(abs((CC-CC_) ./ CC_),[],'all'));
        maxDev(4) = max(maxDev(4),max(abs((MSE-MSE_) ./ MSE_),[],'all'));
    end
end

maxDev
%
%