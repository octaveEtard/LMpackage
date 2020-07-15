% Basic example demonstrating how to use the LM package to fit a linear
% forward model reconstructing the EEG data based on a stimulus.
%
% The data used in this example is synthetic (see: 'generate_example_data')
% giving us access to the ground truth

% folder in which the data is stored ; make sure to run
% 'generate_example_data' before hand, and edit folder name accordingly
baseDataFolder = 'testing_data_2020_07_08_15_21_11';

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


% Time region in which to derive the TRF. Time lag is understood as lag of
% predictor (here stimulus) with respect to predicted data (here EEG).
% Hence, here negative time lags correspond to the causal part of the TRF
% (stimulus preceding response).
minLagT = -500e-3;
maxLagT = 250e-3;


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
                sprintf('sub_%i_cond_%i_part_%i.set',iSub,iCond,iPart)};
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
% size(XtX) = [nFeatures*nLags,nFeature*nLags,nStimPerFile,size(stimOpt)]
% size(Xty) = [nFeatures*nLags,nChan,nStimPerFile,size(EEGOpt)]
%
% that is, in this case:
%
% size(XtX) = [nFeatures*nLags,nFeature*nLags,nStimPerFile,nCond,nParts]
% size(Xty) = [nFeatures*nLags,nChan,nStimPerFile,nCond,nParts,nSub]
%
% This would correspond to fitting model for each subject and each
% stimulus, if the data was used as is.
%
% if opt.sumSub == true, Xty will be summed over subjects: this corresponds
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
opt.sumFrom = 0; % sum over parts but not condition

% Here XtX & Xty will be of size:
% size(XtX) = [nFeatures*nLags,nFeature*nLags,nCond]
% size(Xty) = [nFeatures*nLags,nChan,nCond]
%
% i.e. we fit one model for each condition and all subjects (models fitted
% using the data from all subjects and all stimuli belonging to the same
% condition).

% false: the predictor data (here stimulus) will be zeros padded at its
% edges. true: no padding.
opt.unpad.do = false;
% removing means = fitting models without offsets
opt.removeMean = true;

% getting XtX and Xty required for model fitting here:
[XtX,Xty] = LM_crossMatrices(stimOpt,EEGopt,opt,'forward');

% options to fit he model
trainOpt = struct();
trainOpt.printOut = false;
trainOpt.accumulate = true; % input are XtX & Xty, and not X and Y
trainOpt.method.name = 'ridge-eig-XtX'; % ridge regression
 % normalised regularisation parameter
trainOpt.method.lambda = 10.^(-6:0.5:6);
trainOpt.method.removeEig.type = 'tol';
trainOpt.method.removeEig.tol = 0; % default tolerance

% same stimuli for all subjects
XtX = nSub * XtX;

nLambda = numel(trainOpt.method.lambda);
nLags = opt.maxLag - opt.minLag + 1;

% fit the model independently for each condition, and store the obtained
% coefficients
coeffs = nan(nLags,nChan,nLambda,nCond);
for iCond = 1:nCond
    model = LM_fitLinearModel(XtX(:,:,iCond),Xty(:,:,iCond),trainOpt);
    coeffs(:,:,:,iCond) = model.coeffs; 
end


%% Plotting the obtained TRF against ground truth TRF for comparison
% time vector associated with the coefficients;
tms = 1e3 * LM_getTime(opt,Fs,'forward');

% regularisation coefficient for which to plot the coeffcients
lambda0 = 1e-6;
[~,iLambda0] = min(abs(trainOpt.method.lambda-lambda0));

figure;
for iCond = 1:nCond
    
    % loading ground truth TRF
    trueIR = load(fullfile(baseDataFolder,'IR',sprintf('IR_cond_%i.mat',iCond)));
    tIR = 1e3 * trueIR.tIR;
    trueIR = trueIR.impResponse;
    
    ax = subplot(nCond,1,iCond); hold on;
    c = mean(coeffs(:,:,iLambda0,iCond),2); % mean over channels
    plot(tIR,trueIR);
    plot(tms,c);
    
    ax.XAxis.Label.String = 'Time (ms)';
    ax.YAxis.Label.String = 'Ampltitude (a.u.)';
    
    legend(ax,{'True TRF','Computed TRF'},'Box','off','Location','northwest');
end