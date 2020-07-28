% Basic example demonstrating how to use the LM package to fit a linear
% forward model reconstructing the EEG data based on a stimulus.
%

allSID = {'YH00','YH01','YH02'};%,'YH03','YH04','YH06','YH07','YH08','YH09','YH10','YH11','YH14','YH15','YH16','YH17','YH18','YH19','YH20'};

parts = 1:4;
condition = 'clean';
Fs = 100;

nChan = 64;
procEEG = 'BP-1-12-INTP-AVR';

typeEnv = 'rectified';
procEnv = 'LP-12';


% Time region in which to derive the TRF. Time lag is understood as lag of
% predictor (here stimulus) with respect to predicted data (here EEG).
% Hence, here negative time lags correspond to the causal part of the TRF
% (stimulus preceding response).
minLagT = -1500e-3;
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

% load channel location
chanLocs = LM_example_loadChanLocs();
% define a channel order to be used for all data
chanOrder = {chanLocs(:).labels};

for iPart = 1:nParts
    envFileName = sprintf('env_Fs-%i-%s-%s_%s_%i.mat',Fs,procEnv,typeEnv,condition,iPart);
    stimOpt{iPart} = fullfile(pwd(),'envelopes',envFileName);
    
    for iSub = 1:nSub
        EEGFolder = fullfile(pwd(),'EEGdata',allSID{iSub});
        EEGFileName = sprintf('%s-Fs-%i-%s_%s_%i.set',procEEG,Fs,allSID{iSub},condition,iPart);
        
        EEGopt{iPart,iSub} = {EEGFolder,EEGFileName,chanOrder};
    end
end

% options passed to the call to get the appropriate matrices to fit the
% linear model
opt = struct();
opt.nStimPerFile = 1;
% These are loading function taking one element of stimOpt and EEGopt
% respectively as input, and loading stimulus / EEG data.
opt.getStimulus = @LM_loadFeature_octave;
% This function should return as 1st output a [nPnts x nChan] data matrix,
% and as 2nd outut a vector of indices (size nStimPerFile x 1) indicating
% where each stimulus begins in the data. These indices should be sorted in
% the same order as the stimuli returned by opt.getStimulus.
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
opt.sumFrom = 0; % get all dimenions out

% false: the predictor data (here stimulus) will be zeros padded at its
% edges. true: no padding.
opt.unpad.do = false;
% removing means = fitting models without offsets
opt.removeMean = true;

% convert to samples
opt.nPntsPerf = ceil(tWinPerf*Fs)+1;

% getting XtX and Xty required for model fitting here:
% Here XtX & Xty will be of size:
% size(XtX) = [nLags,nLags,nParts]
% size(Xty) = [nLags,nChan,nParts,nSub]
[XtX,Xty] = LM_crossMatrices(stimOpt,EEGopt,opt,'forward');

% options to fit the model
trainOpt = struct();
trainOpt.method.name = 'ridge-eig-XtX'; % use ridge regression
% regularisation coefficients for which we'll fit the model
trainOpt.method.lambda = 10.^(-6:0.1:6);
trainOpt.method.normaliseLambda = true;
trainOpt.accumulate = true; % the input is XtX & Xty, and not X & y

nLambda = numel(trainOpt.method.lambda);
nLags = opt.maxLag - opt.minLag + 1;
nPerfSize = numel(tWinPerf);

% We will use a leave-one-part-out and leave-one-subject-out
% cross-validation procedure. A model will be fitted over all the data bar
% one part and for all subject bar one. The excluded part for the exluded
% subject will be used as testing data.

% for testing
opt.predBatchSize = nLambda;

CC = cell(nPerfSize,nParts,nSub);
MSE = cell(nPerfSize,nParts,nSub);

for iTestPart = 1:nParts
    
    iTrainParts = [1:(iTestPart-1),(iTestPart+1):nParts];
    XtX_train = sum(XtX(:,:,iTrainParts),3);
    
    % For forward models where all subjects had the same stimuli, a generic
    % model trained over all subjects & stimuli is the same as the average
    % of all models trained for all stimuli and each subject independently
    %
    % We'll take advantage of this here, and train one model for all the
    % subjects, by expanding Xty along the "channel" dimenions (since
    % channels are fitted independently from each other). This avoid
    % having to invert XtX_train for each subject, but requires more memory.
    % Alternatively, see below for the straightforward version.
    Xty_train = reshape( sum(Xty(:,:,iTrainParts,:),3), nLags,nChan*nSub );
    model = LM_fitLinearModel(XtX_train,Xty_train,trainOpt);
    model = reshape(model.coeffs,[nLags,nChan,nSub,nLambda]);
    
    for iTestSub = 1:nSub
        
        iTrainSub = [1:(iTestSub-1),(iTestSub+1):nSub];
        model_train = squeeze(mean(model(:,:,iTrainSub,:),3));
        
        % --- alternatively to the above, do:
        % Xty_train = sum(Xty(:,:,iTrainParts,iTrainSub),[3,4]);
        % model_train = LM_fitLinearModel((nSub-1)*XtX_train,Xty_train,trainOpt);
        % model_train = model_train.coeffs;
        % ---
        
        stim_test = stimOpt(iTestPart);
        EEG_test = EEGopt(iTestPart,iTestSub);
        
        [ CC(:,iTestPart,iTestSub),...
            MSE(:,iTestPart,iTestSub) ] = LM_testModel(model_train,stim_test,EEG_test,opt,'forward');
    end
end

% finally fit one model over all the data
XtX = nSub * sum(XtX,3);
Xty = sum(Xty,[3,4]);

% CC / MSE can be used to choose an appropriate regularisation coefficient,
% and the model obtained for this regularisation coeffcient can be used
% e.g. on held out data.
model = LM_fitLinearModel(XtX,Xty,trainOpt);

%%
% looking at the data using 10s slices
dur0 = 10;
iDur0 = find(tWinPerf == dur0,1);

CC0 = vertcat(CC{iDur0,:});
MSE0 = vertcat(MSE{iDur0,:});

mCC = squeeze(mean(CC0,1))';
mMSE = squeeze(mean(MSE0,1))';

% regularisation curve for each channel
figure;
ax = subplot(1,2,1);
plot(trainOpt.method.lambda,mCC);
ax.XAxis.Scale = 'log';
ax.XAxis.Label.String = '\lambda_n';
ax.YAxis.Label.String = 'Correlation coefficient';

ax = subplot(1,2,2);
plot(trainOpt.method.lambda,mMSE);
ax.XAxis.Scale = 'log';
ax.XAxis.Label.String = '\lambda_n';
ax.YAxis.Label.String = 'MSE';


%% Plotting the obtained TRF against ground truth TRF for comparison
% time vector associated with the coefficients;
tms = 1e3 * LM_getTime(opt,Fs,'forward');

% regularisation coefficient for which to plot the coeffcients
lambda0 = 1e0;
[~,iLambda0] = min(abs(trainOpt.method.lambda-lambda0));

figure;

ax = axes(); hold on;
c = model.coeffs(:,:,iLambda0);
plot(tms,c);

ax.XAxis.Label.String = 'Time (ms)';
ax.YAxis.Label.String = 'Ampltitude (a.u.)';
ax.Title.String = 'TRF';

%% Plotting topography
t0 = 150; % ms 
[~,it0] = min(abs(tms-t0));

figure;
topoplot(model.coeffs(it0,:,iLambda0),chanLocs);
title(sprintf('TRF topography at t = %i ms',t0))
