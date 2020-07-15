% Basic example demonstrating how to use the LM package to fit a linear
% backward model reconstructing the stimulus based on EEG data.
%

allSID = {'YH00','YH01','YH02','YH03','YH04','YH06','YH07','YH08','YH09','YH10','YH11','YH14','YH15','YH16','YH17','YH18','YH19','YH20'};

parts = 1:4;
condition = 'clean';
Fs = 100;

nChan = 64;
procEEG = 'BP-1-12-INTP-AVR';

typeEnv = 'rectified';
procEnv = 'LP-12';

% Time region in which to derive the decoder. Time lag is understood as lag of
% predictor (here EEG) with respect to predicted data (here stimulus).
% Hence, here positive time lags correspond to the causal part of the
% decoder (response after stimulus).
minLagT = -100e-3;
maxLagT = 500e-3;

% estimate performance (CC & MSE) on windows of this duration
% (negative = use all available data)
tWinPerf = [5,10,25,-1]; % in seconds

% stimOpt and EEGopt are arbitray multi-dimensional variables (e.g.
% matrices, cell, structures) containing the required information to load
% stimulus and EEG data respectively. Each element of stimOpt and EEGopt
% will be passed to user-defined loading function (see below) that will
% load the data accordingly. Hence arbitray parameters can be passed to the
% loading functions.
%
% In this example, we will simply store the path of the files to load in
% cell arrays.
nParts = numel(parts);
nSub = numel(allSID);
stimOpt = cell(nParts,1);
% EEGopt have to be of size: [size(stimOpt),nSub], i.e. each stimulus file
% corresponds to 1 EEG recording per subject.
EEGopt = cell(nParts,nSub);

for iPart = 1:nParts
    envFileName = sprintf('env_Fs-%i-%s-%s_%s_%i.mat',Fs,procEnv,typeEnv,condition,iPart);
    stimOpt{iPart} = fullfile(pwd(),'envelopes',envFileName);
    
    for iSub = 1:nSub
        EEGFolder = fullfile(pwd(),'EEGdata',allSID{iSub});
        EEGFileName = sprintf('%s-Fs-%i-%s_%s_%i.set',procEEG,Fs,allSID{iSub},condition,iPart);
        
        EEGopt{iPart,iSub} = {EEGFolder,EEGFileName};
    end
end

% options passed to the call to get the appropriate matrices to fit the
% linear model
opt = struct();
opt.nStimPerFile = 1;
% These are loading function taking one element of stimOpt and EEGopt
% respectively as input, and loading stimulus / EEG data.
opt.getStimulus = @LM_loadFeature_octave;
opt.getResponse = @LM_loadEEG_octave;

% nb of features describing each stimulus
opt.nFeatures = 1;
% nb of channels in the EEG data
opt.nChan = nChan;

% converting lags for time to indices
opt.minLag = floor(minLagT * Fs);
opt.maxLag = ceil(maxLagT * Fs);

opt.sumSub = false;
opt.sumStim = false; % does not matter here
opt.sumFrom = 1; % sum over parts

% false: the predictor data (here stimulus) will be zeros padded at its
% edges. true: no padding.
opt.unpad.do = false;
% removing means = fitting models without offsets
opt.removeMean = true;

% convert to samples
opt.nPntsPerf = ceil(tWinPerf*Fs)+1;


nLags = opt.maxLag - opt.minLag + 1;

% options to fit he model
trainOpt = struct();
trainOpt.printOut = false;
trainOpt.accumulate = true; % input are XtX & Xty, and not X and Y
trainOpt.method.name = 'ridge-eig-XtX'; % ridge regression
% normalised regularisation parameter
trainOpt.method.lambda = 10.^(-6:0.5:6);
trainOpt.method.removeEig.type = 'tol';

nLambda = numel(trainOpt.method.lambda);
nPerfSize = numel(tWinPerf);

% We will fit subject specific models using a leave-one-part-out
% cross-validation procedure. For each subject, a model will be fitted over
% all the data bar one part. The excluded part subject will be used as
% testing data.

CC = cell(nPerfSize,nParts,nSub);
MSE = cell(nPerfSize,nParts,nSub);

for iTestPart = 1:nParts
    
    iTrainParts = [1:(iTestPart-1),(iTestPart+1):nParts];
    
    for iSub = 1:nSub

        % model fitted using only training parts for iSub
        [XtX_train,Xty_train] = LM_crossMatrices(stimOpt(iTrainParts),EEGopt(iTrainParts,iSub),opt,'backward');

        model_train = LM_fitLinearModel(XtX_train,Xty_train,trainOpt);
        model_train = model_train.coeffs;
        
        % testing on the remaining part
        stim_test = stimOpt(iTestPart);
        EEG_test = EEGopt(iTestPart,iSub);
        
 
        [ CC(:,iTestPart,iSub),...
            MSE(:,iTestPart,iSub)] = LM_testModel(model_train,stim_test,EEG_test,opt,'backward');
    end
end


%%
% looking at the data using 10s slices
dur0 = 10;
iDur0 = find(tWinPerf == dur0,1);

CC0 = vertcat(CC{iDur0,:});
nWin = size(CC0,1) / nSub;
CC0 =  reshape(CC0,[nWin,nSub,nLambda]);
mCC = squeeze(mean(CC0,1))';

% regularisation curve for each subject
figure;
ax = axes();
plot(trainOpt.method.lambda ,mCC);
ax.XAxis.Scale = 'log';
ax.XAxis.Label.String = '\lambda_n';
ax.YAxis.Label.String = 'Correlation coefficient';
ax.Title.String = 'Regularisation curve for each subject';

% best CC for each subject
[maxCC,iMax] = max(mCC,[],1);
stdCC = arrayfun(@(iSub,iMax) std(CC0(:,iSub,iMax),[],1),1:nSub,iMax);

[maxCC,iSort] = sort(maxCC,'ascend');
stdCC = stdCC(iSort);

figure;
ax = axes();
errorbar(1:nSub,maxCC,stdCC,'ko');
ax.XAxis.Label.String = 'Subject #';
ax.YAxis.Label.String = 'Correlation coefficient';
ax.Title.String = 'Sorted best correlation coefficients';
ax.XAxis.Limits = [0,nSub+1];
% The decoder obtained at e.g. the best CC for each subject could then be
% used on some held out data.
%
%