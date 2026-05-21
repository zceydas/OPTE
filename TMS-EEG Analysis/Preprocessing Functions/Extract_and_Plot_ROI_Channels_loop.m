clear;
close all;
clc;

%% Define subjects, sessions, and baseline parameters
usr_list = {'007', '011', '012', '013', '017', '019'};
%, {'005', '007', '011', '012', '013', '017', '019'};
sessions = {'Baseline TMSEEG', 'Dosing Session', '1-week follow-up', '2-week follow-up', '1-month follow-up'};

%% Master Loop Across All Subjects
for u = 1:length(usr_list)
    usr_num = usr_list{u};
    
    %% Master Loop Across All Sessions
    for s = 1:length(sessions)
        sess_name = sessions{s};
        
        % Dynamic Epoch Handling: Define target folders based on current session
        if strcmp(sess_name, 'Dosing Session')
            epoch_list = {'Epoch1', 'Epoch2', 'Epoch3', 'Epoch4'};
        else
            epoch_list = {'Epoch0'};
        end
        
        %% Nested Loop Across Epochs Specific to the Active Session
        for e = 1:length(epoch_list)
            epoch_folder = epoch_list{e};
            
            fprintf('\n======================================================================\n');
            fprintf('PROCESSING: Subject %s | Session: %s | Folder: %s\n', usr_num, sess_name, epoch_folder);
            fprintf('======================================================================\n');
            
            % Build dynamic data path
            dataPath = ['/Volumes/T7/OPTE/' usr_num '/Preprocessed TMS/' sess_name '/' epoch_folder];
            
            % Verify directory existence before proceeding
            if ~isfolder(dataPath)
                fprintf('Warning: Path does not exist. Skipping: %s\n', dataPath);
                continue;
            end
            
            % Identify preprocessed data files ending in .vhdr.mat
            dataFiles = dir(fullfile(dataPath, '*.vhdr.mat'));
            
            % Filter out Apple hidden file artifacts (starting with '._')
            validFileIdx = [];
            for f = 1:length(dataFiles)
                if ~startsWith(dataFiles(f).name, '._')
                    validFileIdx = [validFileIdx, f];
                end
            end
            
            if isempty(validFileIdx)
                fprintf('Warning: No valid .vhdr.mat files found in %s. Skipping.\n', dataPath);
                continue;
            end
            
            % Load the first valid dataset found
            dataFilePath = fullfile(dataPath, dataFiles(validFileIdx(1)).name);
            fprintf('Loading file: %s\n', dataFiles(validFileIdx(1)).name);
            
            try
                load(dataFilePath);
            catch ME
                fprintf('Error loading %s. Moving to next set. Message: %s\n', dataFilePath, ME.message);
                continue;
            end
            
            %% Set up output directory
            outputDir = [dataPath '/Amplitude Analysis/'];
            if ~isfolder(outputDir)
                mkdir(outputDir);
                fprintf('Created directory: %s\n', outputDir);
            end
            
            %% Define multiple ROI configurations (LARGE ROIs)
            roiConfigs = struct();
            roiConfigs(1).name = 'LeftMotor';   roiConfigs(1).elecs = {'C1','C3','C5','Cz'};
            roiConfigs(2).name = 'RightMotor';  roiConfigs(2).elecs = {'C1','C2'};
            roiConfigs(3).name = 'LeftFPN';     roiConfigs(3).elecs = {'F1','F3','F5','FC1','FC3','FC5'};
            roiConfigs(4).name = 'RightFPN';    roiConfigs(4).elecs = {'F1','F2','F3','F4'};
            roiConfigs(5).name = 'BothFPN';     roiConfigs(5).elecs = {'P1','P3'};
            roiConfigs(6).name = 'LeftDMN';     roiConfigs(6).elecs = {'P1','P3','P5','CP1','CP3','CP5'};
            roiConfigs(7).name = 'RightDMN';    roiConfigs(7).elecs = {'P2'};
            
            %% Extract all ROIs using TESA function
            fprintf('Extracting ROIs...\n');
            for r = 1:length(roiConfigs)
                EEG = pop_tesa_tepextract(EEG, 'ROI', 'elecs', roiConfigs(r).elecs, 'tepName', roiConfigs(r).name);
            end
            
            %% Perform region-specific peak analysis
            fprintf('Performing region-specific peak analysis...\n');
            
            % DLPFC/FPN peak analysis
            EEG = pop_tesa_peakanalysis(EEG, 'ROI', 'positive', [60], [45 75], 'method', 'largest', 'samples', 5, 'tepName', 'P60_analysis');

            % 2. Detect N100 (Negative peak, targeting 100ms, window 80-140ms)
            EEG = pop_tesa_peakanalysis(EEG, 'ROI', 'negative', [100], [80 140], 'method', 'largest', 'samples', 5, 'tepName', 'N100_analysis');
          
            
            %% Output peak analysis in tables and plots
            fprintf('Exporting peak analysis outputs...\n');
            
            peakOutputs = struct();
            for r = 1:length(roiConfigs)
                roiName = roiConfigs(r).name;
                peakOutput = pop_tesa_peakoutput(EEG, 'tepName', roiName, 'calcType', 'amplitude', 'winType', 'individual', 'averageWin', [], 'fixedPeak', [], 'tablePlot', 'on');
                peakOutputs(r).name = roiName;
                peakOutputs(r).data = peakOutput;
                
                % Explicitly label file identifiers by Subject and Epoch
                tableFileName = [outputDir 'PeakOutput_' roiName '_Subject' usr_num '.xlsx'];
                try
                    if istable(peakOutput)
                        writetable(peakOutput, tableFileName);
                    elseif isstruct(peakOutput)
                        if isfield(peakOutput, 'table') && istable(peakOutput.table)
                            writetable(peakOutput.table, tableFileName);
                        else
                            writetable(struct2table(peakOutput), tableFileName);
                        end
                    elseif ismatrix(peakOutput)
                        writetable(array2table(peakOutput), tableFileName);
                    else
                        csvFileName = [outputDir 'PeakOutput_' roiName '_Subject' usr_num '.csv'];
                        writetable(table(peakOutput), csvFileName);
                    end
                catch ME
                    fprintf('  Warning: Could not save peak output for %s. Error: %s\n', roiName, ME.message);
                end
            end
            
            %% Calculate peak-to-peak amplitudes for each ROI
            fprintf('Calculating peak-to-peak amplitudes...\n');
            
            peakToPeakConfig = struct();
            peakToPeakConfig(1).name = 'LeftMotor';  peakToPeakConfig(1).channels = {'C1','C3','C5','Cz'};          peakToPeakConfig(1).timeWindow = [10 40];
            peakToPeakConfig(2).name = 'RightMotor'; peakToPeakConfig(2).channels = {'C1','C2'};                    peakToPeakConfig(2).timeWindow = [10 40];
            peakToPeakConfig(3).name = 'LeftFPN';    peakToPeakConfig(3).channels = {'F1','F3','F5','FC1','FC3','FC5'}; peakToPeakConfig(3).timeWindow = [20 60];
            peakToPeakConfig(4).name = 'RightFPN';   peakToPeakConfig(4).channels = {'F1','F2','F3','F4'};          peakToPeakConfig(4).timeWindow = [20 60];
            peakToPeakConfig(5).name = 'BothFPN';    peakToPeakConfig(5).channels = {'P1','P3'};                    peakToPeakConfig(5).timeWindow = [20 60];
            peakToPeakConfig(6).name = 'LeftDMN';    peakToPeakConfig(6).channels = {'P1','P3','P5','CP1','CP3','CP5'}; peakToPeakConfig(6).timeWindow = [20 60];
            peakToPeakConfig(7).name = 'RightDMN';   peakToPeakConfig(7).channels = {'P2'};                         peakToPeakConfig(7).timeWindow = [20 60];
            
            peakToPeakTable = table();
            for p = 1:length(peakToPeakConfig)
                roiName = peakToPeakConfig(p).name;
                channels = peakToPeakConfig(p).channels;
                timeWindow = peakToPeakConfig(p).timeWindow;
                
                if exist('jg_peakToPeakDB', 'file') == 2
                    peakToPeakDB = jg_peakToPeakDB(EEG, channels, timeWindow);
                else
                    peakToPeakDB = NaN;
                end
                
                newRow = table({roiName}, peakToPeakDB, 'VariableNames', {'ROI', 'PeakToPeak (µV)'});
                peakToPeakTable = [peakToPeakTable; newRow];
            end
            
            %% Create summary statistics results table
            resultsTable = table();
            for r = 1:length(roiConfigs)
                roiName = roiConfigs(r).name;
                
                if isfield(EEG.ROI, roiName)
                    roi = EEG.ROI.(roiName);
                    
                    if isfield(roi, 'tseries')
                        tseries = roi.tseries;
                        gMean = mean(abs(tseries));
                        peak = max(abs(tseries));
                        [~, peakIdx] = max(abs(tseries));
                        peakLatency = roi.time(peakIdx);
                        auc = trapz(abs(tseries));
                    else
                        gMean = NaN; peak = NaN; peakLatency = NaN; auc = NaN;
                    end
                    
                    newRow = table({roiName}, gMean, peak, peakLatency, auc, ...
                        'VariableNames', {'ROI', 'GlobalMean (µV)', 'PeakAmplitude (µV)', 'PeakLatency (ms)', 'AUC (µV·ms)'});
                    resultsTable = [resultsTable; newRow];
                end
            end
            
            %% Plot TEPs for all ROIs
            fprintf('Generating TEP plots...\n');
            for roiIdx = 1:length(roiConfigs)
                roiName = roiConfigs(roiIdx).name;
                
                if ~isfield(EEG.ROI, roiName), continue; end
                
                roiChanNames = EEG.ROI.(roiName).chans;
                times = EEG.ROI.(roiName).time;
                
                if iscell(roiChanNames)
                    roiChannelIndices = [];
                    for i = 1:length(roiChanNames)
                        idx = find(strcmp({EEG.chanlocs.labels}, roiChanNames{i}));
                        roiChannelIndices = [roiChannelIndices, idx];
                    end
                else
                    roiChannelIndices = roiChanNames;
                end
                
                numChannels = length(roiChannelIndices);
                roiChannels = roiChanNames;
                if numChannels == 0, continue; end
                
                roiData = EEG.data(roiChannelIndices, :, :);
                
                if numChannels <= 2, numRows = 1; numCols = numChannels;
                elseif numChannels <= 4, numRows = 2; numCols = 2;
                elseif numChannels <= 6, numRows = 2; numCols = 3;
                elseif numChannels <= 8, numRows = 2; numCols = 4;
                else, numRows = ceil(numChannels / 4); numCols = 4;
                end
                
                % 'Visible','off' ensures processing runs cleanly in the background
                figH = figure('Position', [100 100 2000 1400], 'Visible', 'off'); 
                set(figH, 'PaperPositionMode', 'auto');
                
                for i = 1:numChannels
                    subplot(numRows, numCols, i);
                    
                    channelData = roiData(i, :, :);
                    channelData = permute(channelData, [2, 3, 1]);
                    meanTEP = mean(channelData, 2);
                    
                    numEpochs = size(channelData, 2);
                    if numEpochs > 1
                        steTEP = std(channelData, [], 2) / sqrt(numEpochs);
                    else
                        steTEP = zeros(size(meanTEP));
                    end
                    
                    hold on;
                    fill([times, fliplr(times)], [meanTEP' + steTEP', fliplr(meanTEP' - steTEP')], ...
                        'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                    plot(times, meanTEP, 'b', 'LineWidth', 2);
                    plot([times(1) times(end)], [0 0], 'k--', 'LineWidth', 0.5);
                    xline(0, 'r--', 'LineWidth', 0.5);
                    
                    componentTimes = [5, 8, 15, 30, 45, 60, 100];
                    for j = 1:length(componentTimes)
                        xline(componentTimes(j), 'k:', 'LineWidth', 1.0);
                    end
                    
                    xlabel('Time (ms)', 'FontSize', 16);
                    ylabel('Amplitude (µV)', 'FontSize', 16);
                    title(roiChannels{i}, 'FontSize', 12);
                    xlim([-100 350]);
                    grid on;
                    hold off;
                end
                
                sgtitle(sprintf('TMS-Evoked Potentials (TEPs) - %s ROI (%s - %s)', roiName, sess_name, epoch_folder), 'FontSize', 16, 'FontWeight', 'bold');
                
                figFileName = [outputDir 'TEP_' roiName '_Subject' usr_num '.fig'];
                pngFileName = [outputDir 'TEP_' roiName '_Subject' usr_num '.png'];
                savefig(figH, figFileName);
                print(figH, pngFileName, '-dpng', '-r300');
                close(figH);
            end
            
            %% Save Final Session & Epoch Summaries
            summaryFileName = [outputDir 'Summary_Metrics_Subject' usr_num '.xlsx'];
            p2pFileName = [outputDir 'Summary_PeakToPeak_Subject' usr_num '.xlsx'];
            
            try
                writetable(resultsTable, summaryFileName);
                writetable(peakToPeakTable, p2pFileName);
                fprintf('Successfully saved Master Summaries for Subject %s [%s - %s]\n', usr_num, sess_name, epoch_folder);
            catch ME
                fprintf('  Warning: Could not save summary tables. Error: %s\n', ME.message);
            end

            close all;
            
        end % End of Epoch Loop
    end % End of Session Loop
end % End of Subject Loop

fprintf('\n*** Completed batch process loop for all subjects, sessions, and multi-epoch pipelines. ***\n');
