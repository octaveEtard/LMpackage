nSub = 3;
nCond = 2;
nParts = 2;
nStimPerFile = 4;

nChan = 5;
noiseAmp = 5;

Fs = 100;
durStim = 3 * 60;

% will timestamp the data generated
timeStamp = datestr(now,'yyyy_mm_dd_HH_MM_SS');
saveFolder = sprintf('testing_data_%s',timeStamp);

saveFolderEEG = fullfile(saveFolder,'EEG');
if ~exist(saveFolderEEG,'dir')
    mkdir(saveFolderEEG);
end

nEdges = ceil(5 * Fs);

for iCond = 1:nCond
    
    [impResponse,tIR] = LM_testing_makeIR(iCond,Fs);
    
    % --- saving
    d = struct();
    d.impResponse = impResponse;
    d.Fs = Fs;
    d.tIR = tIR;
    
    saveName = sprintf('IR_cond_%i',iCond);
    proper_save(d,saveName,fullfile(saveFolder,'IR'));
    % ---
    
    for iPart = 1:nParts
        
        features = cell(nStimPerFile,1);
        
        for iStim = 1:nStimPerFile
            
            % making the stimuli all a different duration
            durStim_ = durStim + randi(10,1,1);
            features{iStim} = LM_testing_makeFeature_continuous(durStim_,Fs,0.5);
            
            % --- saving
            d = struct();
            d.feature = features{iStim};
            d.Fs = Fs;
            
            saveName = sprintf('feature_cond_%i_part_%i_stim_%i',iCond,iPart,iStim);
            proper_save(d,saveName,fullfile(saveFolder,'features'));
            % ---
        end
        
        for iSub = 1:nSub
            
            stimOrder = randperm(nStimPerFile);
            iB = nan(nStimPerFile,1);
            
            r = cell(nStimPerFile,1);
            
            for iiStim = 1:nStimPerFile
                [r{iiStim},iB(iiStim)] = LM_testing_makeResponse(features{stimOrder(iiStim)},impResponse,nEdges,nChan,noiseAmp);
            end
            
            n = cellfun(@(m) size(m,1), r);
            iB = iB + cumsum([0;n(1:end-1)]);
            
            r = vertcat(r{:});
            
            %%
            EEG = eeg_emptyset;
            EEG.data = r';
            EEG.trials = 1;
            EEG.nbchan = size(EEG.data, 1);
            EEG.pnts = size(EEG.data, 2);
            EEG.event = struct();
            EEG.srate = Fs;
            
            for iiStim = 1:nStimPerFile
                EEG.event(iiStim).latency = iB(iiStim);
                EEG.event(iiStim).duration = 1/Fs;
                EEG.event(iiStim).type = 'Stim begin';
                EEG.event(iiStim).code = stimOrder(iiStim);
                EEG.event(iiStim).urevent = iiStim;
                
            end
            EEG.setname = sprintf('sub_%i_cond_%i_part_%i',iSub,iCond,iPart);
            EEG = eeg_checkset(EEG);
            
            pop_saveset(EEG,'filename',EEG.setname,'filepath',saveFolderEEG);
        end
    end
end
