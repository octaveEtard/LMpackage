% Basic example demonstrating how to use the LM package to fit a linear
% backward model reconstructing the stimulus based on EEG data.
% Author: Octave Etard
%
% The data used in this example is synthetic (see: 'generate_example_data')
% giving us access to the ground truth

% folder in which the data is stored ; make sure to run
% 'generate_example_data' before hand, and edit folder name accordingly
baseDataFolder = 'testing_data_someTimestamp';

% Descritpion of the synthetic data:
d = load(fullfile(baseDataFolder,'generating_parameter.mat'));

nChan = d.nChan; % nb of channels in the EEG
nCond = d.nCond; % nb of conditions
nParts = d.nParts; % the data from each condition is divided in nParts files
nStimPerFile = d.nStimPerFile; % each EEG file contains responses to nStimPerFile stimuli

% nb of subjects
nSub = d.nSub;

% The structure of the data (conditions/parts/stimPerFile...) is arbitray,
% but for each subject analysed here, the EEG files have to contain
% responses to the same set of stimuli. The presentation order of the
% stimuli may be different

Fs = d.Fs; % sampling rate


% Time region in which to derive the decoder. Time lag is understood as lag of
% predictor (here EEG) with respect to predicted data (here stimulus).
% Hence, here positive time lags correspond to the causal part of the
% decoder (response after stimulus).
minLagT = -250e-3;
maxLagT = 500e-3;


% stimOpt and EEGopt are arbitray multi-dimensional variables (e.g.
% matrices, cell, structures) containing the required information to load
% stimulus and EEG data respectively. Each element of stimOpt and EEGopt
% will be passed to user-defined loading function (see below) that will
% load the data accordingly. Hence arbitray parameters can be passed to the
% loading functions.
%
% In this example, we will simply store the path of the files to load in
% cell arrays.
stimOpt = cell(nCond,nParts);
% EEGopt have to be of size: [size(stimOpt),nSub], i.e. each stimulus file
% corresponds to 1 EEG recording per subject.
EEGopt = cell(nCond,nParts,nSub);

for iCond = 1:nCond
    for iPart = 1:nParts
        c = cell(nStimPerFile,1);
        for iStim = 1:nStimPerFile
            c{iStim} = fullfile(baseDataFolder,'features',...
                sprintf('feature_cond_%i_part_%i_stim_%i.mat',iCond,iPart,iStim));
        end
        stimOpt{iCond,iPart} = c;
        
        for iSub = 1:nSub
            EEGopt{iCond,iPart,iSub} = {fullfile(baseDataFolder,'EEG'),...
                sprintf('sub_%i_cond_%i_part_%i.mat',iSub,iCond,iPart)};
        end
    end
end

% options passed to the call to get the appropriate matrices to fit the
% linear model
opt = struct();
opt.nStimPerFile = nStimPerFile;
% These are loading function taking one element of stimOpt and EEGopt
% respectively as input, and loading stimulus / EEG data.
opt.getStimulus = @LM_testing_loadFeature;
opt.getResponse = @LM_testing_loadEEG;

% nb of features describing each stimulus
opt.nFeatures = 1;
% nb of channels in the EEG data
opt.nChan = nChan;

% converting lags for time to indices
opt.minLag = floor(minLagT * Fs);
opt.maxLag = ceil(maxLagT * Fs);

% The call will return  matrices XtX and Xty of size:
% size(XtX) = [nChan*nLags,nChan*nLags,nStimPerFile,size(stimOpt),nSub]
% size(Xty) = [nChan*nLags,nFeatures,nStimPerFile,size(EEGOpt)]
%
% that is, in this case:
%
% size(XtX) = [nChan*nLags,nChan*nLags,nStimPerFile,nCond,nParts,nSub]
% size(Xty) = [nChan*nLags,nFeatures,nStimPerFile,nCond,nParts,nSub]
%
% This would correspond to fitting model for each subject and each
% stimulus, if the data was used as is.
%
% if opt.sumSub == true, XtX & Xty will be summed over subjects: this corresponds
% to fitting one model for all subjects, and for each stimulus.
%
% if opt.sumStim == true, XtX & Xty will be summed over their 3rd dimensions
% (stimPerFile). This corresponds to fitting one model for all stimuli in
% one EEG file.
%
% if 1 <= opt.sumFrom, XtX & Xty will be summed over all
% dimensions >= opt.sumFrom in stimOpt. This corresponds to fitting one
% model for all collapsed dimenions. opt.sumFrom == 1 will sum over all
% dimenions, opt.sumFrom <= 0 will not sum over any dimenions.
opt.sumSub = true;
opt.sumStim = true;
opt.sumFrom = 2; % sum over parts but not condition

% Here XtX & Xty will be of size:
% size(XtX) = [nChan*nLags,nFeatures,nCond]
% size(Xty) = [nChan*nLags,nFeatures,nCond]
%
% i.e. we fit one decoder for each condition and all subjects (models fitted
% using the data from all subjects and all stimuli belonging to the same
% condition).

% false: the predictor data (here stimulus) will be zeros padded at its
% edges. true: no padding.
opt.unpad.do = false;
% removing means = fitting models without offsets
opt.removeMean = true;

% getting XtX and Xty required for model fitting here:
[XtX,Xty] = LM_crossMatrices(stimOpt,EEGopt,opt,'backward');

% options to fit the model
trainOpt = struct();
trainOpt.method.name = 'ridge-eig-XtX'; % use ridge regression
% regularisation coefficients for which we'll fit the model
trainOpt.method.lambda = 10.^(-6:0.1:6);
trainOpt.method.normaliseLambda = true;
trainOpt.accumulate = true; % the input is XtX & Xty, and not X & y

nLambda = numel(trainOpt.method.lambda);
nLags = opt.maxLag - opt.minLag + 1;

% fit the model independently for each condition, and store the obtained
% coefficients
coeffs = nan(nLags,nChan,nLambda,nCond);
for iCond = 1:nCond
    model = LM_fitLinearModel(XtX(:,:,iCond),Xty(:,:,iCond),trainOpt);
    coeffs(:,:,:,iCond) = reshape(model.coeffs,[nLags,nChan,nLambda]);
end

% Return a time vector associated with the coefficients, and make sure the
% causal part of the coefficients are at positive latencies.
[tms,coeffs] = LM_getTime(opt,Fs,'backward',coeffs,1);
tms = 1e3 * tms; % in milliseconds


%% Plotting the decoder 
%
% regularisation coefficient for which to plot the coeffcients
lambda0 = 1e-6;
[~,iLambda0] = min(abs(trainOpt.method.lambda-lambda0));

figure;
for iCond = 1:nCond
    
    ax = subplot(nCond,1,iCond); hold on;
    c = mean(coeffs(:,:,iLambda0,iCond),2); % mean over channels
    plot(tms,c);
    
    ax.XAxis.Label.String = 'Time (ms)';
    ax.YAxis.Label.String = 'Ampltitude (a.u.)';
    
    legend(ax,{'Mean decoder coefficients'},'Box','off','Location','northwest');
end



