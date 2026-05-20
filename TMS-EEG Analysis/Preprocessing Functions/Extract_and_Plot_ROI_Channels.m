clear;
close all;
clc;

%% Define subject number
usr_num = '005';

%% Set up output directory (define early so it's available throughout script)
outputDir = ['/Volumes/T7/OPTE/' usr_num '/Analysis/Amplitude_Analysis/'];

% Create directory if it doesn't exist
if ~isfolder(outputDir)
    mkdir(outputDir);
    fprintf('Created directory: %s\n', outputDir);
end


%% Load preprocessed EEG data
% Load the final preprocessed dataset
dataPath = ['/Volumes/T7/OPTE/' usr_num '/Preprocessed TMS/Baseline TMSEEG/Epoch0_2_12/' usr_num '_Baseline TMSEEG_Subject' usr_num '_Baseline_Epoch0_SinglePulse.vhdr_preICARejection.mat'];
load(dataPath);

%% Define multiple ROI configurations (from ExtractROIs_v2.m - LARGE ROIs)
roiConfigs = struct();
roiConfigs(1).name = 'LeftMotor';
roiConfigs(1).elecs = {'C1','C3','C5','Cz'};

roiConfigs(2).name = 'RightMotor';
roiConfigs(2).elecs = {'C1','C2'};

roiConfigs(3).name = 'LeftFPN';
roiConfigs(3).elecs = {'F1','F3','F5','FC1','FC3','FC5'};

roiConfigs(4).name = 'RightFPN';
roiConfigs(4).elecs = {'F1','F2','F3','F4'};

roiConfigs(5).name = 'BothFPN';
roiConfigs(5).elecs = {'P1','P3'};

roiConfigs(6).name = 'LeftDMN';
roiConfigs(6).elecs = {'P1','P3','P5','CP1','CP3','CP5'};

roiConfigs(7).name = 'RightDMN';
roiConfigs(7).elecs = {'P2'};

%roiConfigs(8).name = 'allChans';
%roiConfigs(8).elecs = 'all';

%% Extract all ROIs using TESA function
fprintf('Extracting ROIs...\n');
for r = 1:length(roiConfigs)
    EEG = pop_tesa_tepextract(EEG, 'ROI', 'elecs', roiConfigs(r).elecs, 'tepName', roiConfigs(r).name);
    fprintf('  %s extracted\n', roiConfigs(r).name);
end

%% Perform region-specific peak analysis
fprintf('\nPerforming region-specific peak analysis...\n');

% DLPFC/FPN peak analysis
EEG = pop_tesa_peakanalysis(EEG, 'ROI', 'positive', [30], [20 60], 'method', 'largest', 'samples', 5);
EEG = pop_tesa_peakanalysis(EEG, 'ROI', 'negative', [30], [20 60], 'method', 'largest', 'samples', 5);

% Motor cortex peak analysis
EEG = pop_tesa_peakanalysis(EEG, 'ROI', 'positive', [32 61 180], [26 38; 56 66; 160 200], 'method', 'largest', 'samples', 10);
EEG = pop_tesa_peakanalysis(EEG, 'ROI', 'negative', [15 50 95], [5 25; 45 55; 85 115], 'method', 'largest', 'samples', 10);

% DMN peak analysis
EEG = pop_tesa_peakanalysis(EEG, 'ROI', 'positive', [25 180], [15 35; 170 190], 'method', 'largest', 'samples', 5);
EEG = pop_tesa_peakanalysis(EEG, 'ROI', 'negative', [55 100], [45 65; 90 110], 'method', 'largest', 'samples', 5);

%% Output peak analysis in tables and plots
fprintf('\nExporting peak analysis outputs...\n');

% Process each ROI and save peak output tables
peakOutputs = struct();
for r = 1:length(roiConfigs)
    roiName = roiConfigs(r).name;
    peakOutput = pop_tesa_peakoutput(EEG, 'tepName', roiName, 'calcType', 'amplitude', 'winType', 'individual', 'averageWin', [], 'fixedPeak', [], 'tablePlot', 'on');
    peakOutputs(r).name = roiName;
    peakOutputs(r).data = peakOutput;
    
    % Save peak output - handle different data types
    tableFileName = [outputDir 'PeakOutput_' roiName '_Subject' usr_num '.xlsx'];
    try
        if istable(peakOutput)
            % If it's already a table, save directly
            writetable(peakOutput, tableFileName);
            fprintf('  Saved peak output for %s to: %s\n', roiName, tableFileName);
        elseif isstruct(peakOutput)
            % If it's a struct, convert to table
            peakOutputTable = struct2table(peakOutput);
            writetable(peakOutputTable, tableFileName);
            fprintf('  Saved peak output for %s to: %s (converted from struct)\n', roiName, tableFileName);
        elseif ismatrix(peakOutput)
            % If it's a matrix, create a table from it
            peakOutputTable = array2table(peakOutput);
            writetable(peakOutputTable, tableFileName);
            fprintf('  Saved peak output for %s to: %s (converted from matrix)\n', roiName, tableFileName);
        else
            % For other types, try to save as CSV
            csvFileName = [outputDir 'PeakOutput_' roiName '_Subject' usr_num '.csv'];
            try
                writetable(table(peakOutput), csvFileName);
                fprintf('  Saved peak output for %s to: %s (as CSV)\n', roiName, csvFileName);
            catch
                fprintf('  Warning: Could not save peak output for %s. Unsupported data type.\n', roiName);
            end
        end
    catch ME
        fprintf('  Warning: Could not save peak output for %s. Error: %s\n', roiName, ME.message);
    end
end
fprintf('Peak analysis outputs completed.\n');

%% Calculate peak-to-peak amplitudes for each ROI
fprintf('\nCalculating peak-to-peak amplitudes...\n');

% Define channels and time windows for peak-to-peak analysis by ROI
peakToPeakConfig = struct();
peakToPeakConfig(1).name = 'LeftMotor';
peakToPeakConfig(1).channels = {'C1','C3','C5','Cz'};
peakToPeakConfig(1).timeWindow = [10 40];

peakToPeakConfig(2).name = 'RightMotor';
peakToPeakConfig(2).channels = {'C1','C2'};
peakToPeakConfig(2).timeWindow = [10 40];

peakToPeakConfig(3).name = 'LeftFPN';
peakToPeakConfig(3).channels = {'F1','F3','F5','FC1','FC3','FC5'};
peakToPeakConfig(3).timeWindow = [20 60];

peakToPeakConfig(4).name = 'RightFPN';
peakToPeakConfig(4).channels = {'F1','F2','F3','F4'};
peakToPeakConfig(4).timeWindow = [20 60];

peakToPeakConfig(5).name = 'LeftDMN';
peakToPeakConfig(5).channels = {'P1','P3','P5','CP1','CP3','CP5'};
peakToPeakConfig(5).timeWindow = [20 60];

peakToPeakConfig(6).name = 'RightDMN';
peakToPeakConfig(6).channels = {'P2'};
peakToPeakConfig(6).timeWindow = [20 60];

% Calculate peak-to-peak for each ROI
peakToPeakTable = table();
for p = 1:length(peakToPeakConfig)
    roiName = peakToPeakConfig(p).name;
    channels = peakToPeakConfig(p).channels;
    timeWindow = peakToPeakConfig(p).timeWindow;
    
    if exist('jg_peakToPeakDB', 'file') == 2
        peakToPeakDB = jg_peakToPeakDB(EEG, channels, timeWindow);
        fprintf('  %s: Peak-to-peak = %.2f µV\n', roiName, peakToPeakDB);
    else
        % Fallback: calculate manually if function not available
        fprintf('  Warning: jg_peakToPeakDB not found. Using manual calculation for %s\n', roiName);
        peakToPeakDB = NaN;
    end
    
    % Add to results table
    newRow = table({roiName}, peakToPeakDB, ...
        'VariableNames', {'ROI', 'PeakToPeak (µV)'});
    peakToPeakTable = [peakToPeakTable; newRow];
end

%% Create results table for all ROIs
resultsTable = table();
roiNames = {roiConfigs.name};

for r = 1:length(roiConfigs)
    roiName = roiConfigs(r).name;
    
    if isfield(EEG.ROI, roiName)
        roi = EEG.ROI.(roiName);
        
        % Compute metrics from tseries data
        if isfield(roi, 'tseries')
            tseries = roi.tseries;  % Grand average time series
            
            % Global mean amplitude
            gMean = mean(abs(tseries));
            
            % Peak amplitude (max absolute value)
            peak = max(abs(tseries));
            
            % Peak latency (time of maximum)
            [~, peakIdx] = max(abs(tseries));
            peakLatency = roi.time(peakIdx);
            
            % Area under curve (integral of absolute value)
            auc = trapz(abs(tseries));
        else
            gMean = NaN;
            peak = NaN;
            peakLatency = NaN;
            auc = NaN;
        end
        
        % Add row to results table
        newRow = table({roiName}, gMean, peak, peakLatency, auc, ...
            'VariableNames', {'ROI', 'GlobalMean (µV)', 'PeakAmplitude (µV)', 'PeakLatency (ms)', 'AUC (µV·ms)'});
        resultsTable = [resultsTable; newRow];
    end
end

%% Plot TEPs for all ROIs
fprintf('\nGenerating TEP plots for all ROIs...\n');

for roiIdx = 1:length(roiConfigs)
    roiName = roiConfigs(roiIdx).name;
    
    % Check if ROI was extracted
    if ~isfield(EEG.ROI, roiName)
        fprintf('  Skipping %s - not extracted\n', roiName);
        continue;
    end
    
    %% Get the extracted ROI data
    roiChanNames = EEG.ROI.(roiName).chans;
    times = EEG.ROI.(roiName).time;
    
    %% Convert channel names to indices
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
    
    % Skip if no channels
    if numChannels == 0
        fprintf('  Skipping %s - no channels found\n', roiName);
        continue;
    end
    
    %% Extract individual channel data from original EEG
    roiData = EEG.data(roiChannelIndices, :, :);
    
    %% Create figure based on number of channels
    % Determine subplot layout
    if numChannels <= 2
        numRows = 1; numCols = numChannels;
    elseif numChannels <= 4
        numRows = 2; numCols = 2;
    elseif numChannels <= 6
        numRows = 2; numCols = 3;
    elseif numChannels <= 8
        numRows = 2; numCols = 4;
    else
        numRows = ceil(numChannels / 4); numCols = 4;
    end
    
    figure('Position', [100 100 2000 1400]);
    set(gcf, 'PaperPositionMode', 'auto');
    
    for i = 1:numChannels
        subplot(numRows, numCols, i);
        
        % Get data for this channel (roiData is channels x timepoints x epochs)
        channelData = squeeze(roiData(i, :, :));
        
        % Ensure proper dimensions (timepoints x epochs)
        if size(channelData, 1) < size(channelData, 2)
            channelData = channelData';
        end
        
        % Calculate mean TEP across epochs
        meanTEP = mean(channelData, 2);
        
        % Calculate standard error
        steTEP = std(channelData, [], 2) / sqrt(size(channelData, 2));
        
        % Plot mean TEP with shaded error
        hold on;
        fill([times, fliplr(times)], [meanTEP' + steTEP', fliplr(meanTEP' - steTEP')], ...
            'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(times, meanTEP, 'b', 'LineWidth', 2);
        
        % Add zero line
        plot([times(1) times(end)], [0 0], 'k--', 'LineWidth', 0.5);
        xline(0, 'r--', 'LineWidth', 0.5);
        
        % Add component timing lines (vertical dotted lines)
        componentTimes = [5, 8, 15, 30, 45, 60, 100];
        for j = 1:length(componentTimes)
            xline(componentTimes(j), 'k:', 'LineWidth', 1.0);
        end
        
        % Labels and formatting
        xlabel('Time (ms)', 'FontSize', 16);
        ylabel('Amplitude (µV)', 'FontSize', 16);
        title(roiChannels{i}, 'FontSize', 12);
        xlim([-100 350]);
        grid on;
        set(gca, 'FontSize', 10);
        hold off;
    end
    
    % Overall figure title
    sgtitle(sprintf('TMS-Evoked Potentials (TEPs) - %s ROI', roiName), 'FontSize', 16, 'FontWeight', 'bold');
    
    % Save figure to output directory
    figFileName = [outputDir 'TEP_' roiName '_Subject' usr_num '.fig'];
    pngFileName = [outputDir 'TEP_' roiName '_Subject' usr_num '.png'];
    savefig(figFileName);
    % Save as high-resolution PNG
    print(gcf, pngFileName, '-dpng', '-r300');
    fprintf('  Generated plot for %s (%d channels)\n', roiName, numChannels);
    fprintf('  Saved figures to: %s and %s\n', figFileName, pngFileName);
end

%% Print comprehensive analysis summary
fprintf('\n========== Comprehensive TEP Analysis Summary ==========\n');
fprintf('Subject: 005\n');
fprintf('Session: Baseline TMSEEG\n');

fprintf('\n--- ALL ROIs METRICS ---\n');
disp(resultsTable);

fprintf('\n--- PEAK-TO-PEAK ANALYSIS ---\n');
disp(peakToPeakTable);

fprintf('\n--- INDIVIDUAL ROI CHANNEL STATISTICS ---\n');
for roiIdx = 1:length(roiConfigs)
    roiName = roiConfigs(roiIdx).name;
    
    if ~isfield(EEG.ROI, roiName)
        continue;
    end
    
    roiChanNames = EEG.ROI.(roiName).chans;
    
    if iscell(roiChanNames)
        roiChannelIndices = [];
        for i = 1:length(roiChanNames)
            idx = find(strcmp({EEG.chanlocs.labels}, roiChanNames{i}));
            roiChannelIndices = [roiChannelIndices, idx];
        end
    else
        roiChannelIndices = roiChanNames;
    end
    
    if isempty(roiChannelIndices)
        continue;
    end
    
    roiData = EEG.data(roiChannelIndices, :, :);
    roiChannels = roiChanNames;
    
    fprintf('\n%s ROI:\n', roiName);
    fprintf('  Number of epochs: %d\n', size(roiData, 3));
    for i = 1:length(roiChannelIndices)
        channelData = squeeze(roiData(i, :, :));
        
        if size(channelData, 1) < size(channelData, 2)
            channelData = channelData';
        end
        
        meanAmplitude = mean(channelData, 'all');
        peakAmplitude = max(abs(channelData), [], 'all');
        fprintf('    %s - Mean: %.2f µV, Peak: %.2f µV\n', roiChannels{i}, meanAmplitude, peakAmplitude);
    end
end

fprintf('======================ANALYSIS COMPLETE======================\n\n');
