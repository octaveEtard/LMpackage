% Basic example demonstrating how to use the LM package to fit a linear
% backward model reconstructing the stimulus based on EEG data.
%
idxSID = [2,3]; % subject indices
% subject IDs
allSID = arrayfun(@(i) sprintf('P%02i',i),idxSID,'UniformOutput',false);

%
% In this example each EEG file contains responses to 15 different stories,
% listed below:
stories = {...
    'AUNP01','AUNP02','AUNP03','AUNP04','AUNP05','AUNP06','AUNP07','AUNP08',...
    'BROP01','BROP02','BROP03','FLOP01','FLOP02','FLOP03','FLOP04'};
%
% The matrix below contains the time of stimulus onset for each of the
% stories and each of the subject (matrix of size nSub x nStories).
% Stimuli are sorted in alphabetical order, that is the same order as the
% 'stories' variable above.
alignmentFilePath = fullfile(pwd(),'stimuli','all_story_onsets.mat');

nChan = 64;
Fs = 100;
procEEG = 'BP-1-12-INTP-AVR';

typeEnv = 'rectified';
procEnv = 'BP-1-12';


% Time region in which to derive the decoder. Time lag is understood as lag of
% predictor (here EEG) with respect to predicted data (here stimulus).
% Hence, here positive time lags correspond to the causal part of the
% decoder (response after stimulus).
minLagT = -100e-3;
maxLagT = 500e-3;

% estimate performance (CC & MSE) on windows of this duration
% (negative = use all available data)
tWinPerf = 10; % in seconds

% stimOpt and EEGopt are arbitray multi-dimensional variables (e.g.
% matrices, cell, structures) containing the required information to load
% stimulus and EEG data respectively. Each element of stimOpt and EEGopt
% will be passed to user-defined loading function (see below) that will
% load the data accordingly. Hence arbitray parameters can be passed to the
% loading functions.
%
% In this example, we will simply store the path of the files to load in
% cell arrays.
nStories = numel(stories);
nSub = numel(allSID);

stimOpt = {cell(nStories,1)};
% one EEG file per subjects
EEGopt = cell(nSub,1);


for iStory = 1:nStories
    envFileName = sprintf('%s-Fs-%i-%s-%s.mat',procEnv,Fs,typeEnv,stories{iStory});
    stimOpt{1}{iStory} = fullfile(pwd(),'stimuli',envFileName);
end

for iSub = 1:nSub
    EEGFolder = fullfile(pwd(),'EEGdata',allSID{iSub});
    EEGFileName = sprintf('%s-Fs-%i-%s.set',procEEG,Fs,allSID{iSub});
    
    EEGopt{iSub} = {EEGFolder,EEGFileName,alignmentFilePath,idxSID(iSub),1:nStories};
end

% options passed to the call to get the appropriate matrices to fit the
% linear model
opt = struct();
opt.nStimPerFile = nStories;
% These are loading function taking one element of stimOpt and EEGopt
% respectively as input, and loading stimulus / EEG data.
opt.getStimulus = @LM_loadFeature_hugo;
% This function should return as 1st output a [nPnts x nChan] data matrix,
% and as 2nd outut a vector of indices (size nStimPerFile x 1) indicating
% where each stimulus begins in the data. These indices should be sorted in
% the same order as the stimuli returned by opt.getStimulus.
opt.getResponse = @LM_loadEEG_hugo;

% nb of features describing each stimulus
opt.nFeatures = 1;
% nb of channels in the EEG data
opt.nChan = nChan;

% converting lags for time to indices
opt.minLag = floor(minLagT * Fs);
opt.maxLag = ceil(maxLagT * Fs);

opt.sumSub = false;
opt.sumStim = false; % does not matter here
opt.sumFrom = 0;

% false: the predictor data (here stimulus) will be zeros padded at its
% edges. true: no padding.
opt.unpad.do = false;
% removing means = fitting models without offsets
opt.removeMean = true;

% convert to samples
opt.nPntsPerf = ceil(tWinPerf*Fs)+1;


nLags = opt.maxLag - opt.minLag + 1;

% options to fit the model
trainOpt = struct();
trainOpt.method.name = 'ridge-eig-XtX'; % use ridge regression
% regularisation coefficients for which we'll fit the model
trainOpt.method.lambda = 10.^(-6:0.1:6);
trainOpt.method.normaliseLambda = true;
trainOpt.accumulate = true; % the input is XtX & Xty, and not X & y

nLambda = numel(trainOpt.method.lambda);
nPerfSize = numel(tWinPerf);

% We will fit subject specific models using a leave-one-story-out
% cross-validation procedure. For each subject, a model will be fitted over
% all the data bar one story. The excluded story will be used as
% testing data.

CC = cell(nPerfSize,nStories,nSub);
MSE = cell(nPerfSize,nStories,nSub);

opt_test = opt;
opt_test.nStimPerFile = 1;



for iSub = 1:nSub
    
    % matrices for all data parts for iSub
    [XtX,Xty] = LM_crossMatrices(stimOpt,EEGopt(iSub),opt,'backward');
    
    for iTestStory = 1:nStories
        
        % --- training
        idxTrain = [1:(iTestStory-1),(iTestStory+1):nStories];
        
        XtX_train = sum(XtX(:,:,idxTrain),3);
        Xty_train = sum(Xty(:,:,idxTrain),3);
        
        model = LM_fitLinearModel(XtX_train,Xty_train,trainOpt);
        model = model.coeffs;
        
        % --- testing on the remaining part
        stim_test = {stimOpt{1}(iTestStory)};
        EEG_test = {{EEGopt{iSub}{1:4},iTestStory}};

        [ CC(:,iTestStory,iSub),...
            MSE(:,iTestStory,iSub)] = LM_testModel(model,stim_test,EEG_test,opt_test,'backward');
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