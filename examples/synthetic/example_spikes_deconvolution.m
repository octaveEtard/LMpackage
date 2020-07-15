%% Generate synthetic data as a convolution between an Impulse Response
% (IR) / Temporal Response Function (TRF) & some events (stimulus) and
% then try to recover the IR / TRF from the data & stimulus
% (FORWARD model)
%
tMax = 1 * 60; % data duration, in s
Fs = 100; % sampling rate, in Hz

nPnts = ceil(tMax*Fs) + 1; % number of points
tMax = (nPnts-1)/Fs; % making sure tMax & nPnts match exactly

tEEG = ((1:nPnts)-1)/Fs;


%% Creating vector of events (spikes)
%
tISI = 250e-3; % in s, spacing betwen events
jitter = 50e-3; % in s, event jiter

feature = LM_testing_makeFeature_spikes(tMax,Fs,tISI,jitter);


%% Ground truth response
% Choose between:

% --- simple Gaussian pulse response
type = 1;
%
% OR
% --- bi-modal response
% type = 2;

[impResponse,tResponse,nMaxResponse] = LM_testing_makeIR(type,Fs);


%% Creating the synthetic EEG / response
noiseAmp = 0.5; % noise amplitude (1 = same as max of IR)

[eeg,iB] = LM_testing_makeResponse(feature,impResponse,0,1,0);

% 1 st col: ideal EEG, no noise
% 2 nd col: EEG + noise
eeg = [...
    eeg, LM_testing_makeResponse(feature,impResponse,0,1,noiseAmp)];


%% Deconvolving the stimulus from the EEG data to uncover the neural
% response

% time frame in which we want to derive the neural response
minLag = -nMaxResponse;
maxLag = nMaxResponse;


% forming XtX and Xty:
opt = struct();
opt.removeMean = true;
opt.unpad = LM_laggedDims(size(feature,1),iB,size(eeg,1),minLag,maxLag);
opt.unpad.do = false;

[XtX,xF,mX,Xtop,Xbottom] = LM_laggedXtX(feature,minLag,maxLag,opt);
opt.iB = iB;
opt.nx = size(feature,1);

Xty = LM_laggedXty(xF,eeg,minLag,maxLag,...
    mX,Xtop,Xbottom,...
    false,[],[],[],[],...
    opt);

% options to fit the model
trainOpt = struct();
trainOpt.method.name = 'ridge-eig-XtX';
% regularisation coefficients for which we'll fit the model
trainOpt.method.lambda = 10.^(-6:0.1:6);
trainOpt.method.removeEig.type = 'tol';
trainOpt.method.removeEig.tol = 0; % default tolerance

trainOpt.accumulate = true;
trainOpt.printOut = false;

% fitting the model
model = LM_fitLinearModel(XtX,Xty,trainOpt);


%% Plotting the results
lambda0 = 1e-6; % regularisation coefficient to use;

[~,iLambda0] = min(abs(lambda0 - trainOpt.method.lambda));
lambda0 = trainOpt.method.lambda(iLambda0);

tms = 1e3 * (-maxLag:-minLag) / Fs;
nLambda = numel(trainOpt.method.lambda);
nChan = size(eeg,2);


nLags = maxLag - minLag + 1;
coeffs = reshape(model.coeffs,[nLags,nChan,nLambda]);

t = {'No noise','With noise'};

fig = figure;
for iChan = 1:nChan
    ax = subplot(1,nChan,iChan); hold on;
    plot(1e3*tResponse,impResponse,'k-');
    plot(tms,squeeze(coeffs(:,iChan,iLambda0)));

    ax.Title.String = t{iChan};
    
    legend(ax,{'True IR',sprintf('Recovered IR (log(\\lambda) = 10^{%.1f})',log10(lambda0))},'box','off','Location','northwest');
end
