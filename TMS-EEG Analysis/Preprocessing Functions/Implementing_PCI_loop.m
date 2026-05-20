clear;
close all;
clc;

%% Define subjects, sessions, and parameters
usr_list = {'005', '007', '011', '012', '013', '017', '019'};
sessions = {'Baseline TMSEEG', 'Dosing Session', '1-week follow-up', '2-week follow-up', '1-month follow-up'};

% Base path for your data
baseDataPath = '/Volumes/T7/OPTE/';

% Define PCIst parameters based on Comolatti et al. (2019)
params = struct(...
    'baseline', [-400, -20], ...  % Baseline window (ms)
    'response', [15, 300], ...    % Response window (ms)
    'max_var', 99, ...           % Retains 99% variance
    'min_snr', 1.1, ...          % Minimum SNR threshold
    'k', 1.2, ...                % Default threshold multiplier
    'l', 1, ...                  % Default transition parameter
    'nsteps', 100 ...            % Default integration steps
);

%% Initialize Results Table
% Pre-allocating an empty table to store findings cleanly
results = table();
rowIdx = 1;

%% Nested Loops to iterate through Users, Sessions, and Epochs
for u = 1:length(usr_list)
    current_user = usr_list{u};
    
    for s = 1:length(sessions)
        current_session = sessions{s};
        
        % Determine epochs based on the session name
        if strcmp(current_session, 'Dosing Session')
            epochs_to_run = [1, 2, 3, 4];
        else
            epochs_to_run = 0;
        % You could also use: if s == 2
        end
        
        for e = 1:length(epochs_to_run)
            current_epoch = epochs_to_run(e);
            
            % Construct the directory path dynamically
            dataPath = fullfile(baseDataPath, current_user, 'Preprocessed TMS', current_session, ['Epoch' num2str(current_epoch)]);
            
            % Check if the directory exists before proceeding
            if ~exist(dataPath, 'dir')
                warning('Directory does not exist, skipping: %s', dataPath);
                continue;
            end
            
            % Find the .vhdr.mat files
            dataFiles = dir(fullfile(dataPath, '*.vhdr.mat'));
            if isempty(dataFiles)
                warning('No preprocessed data files found in: %s', dataPath);
                continue;
            end
            
            % Pick the first file and clean macOS hidden file prefix if present
            fileName = dataFiles(1).name;
            if startsWith(fileName, '._')
                fileName = strrep(fileName, '._', '');
            end
            dataFilePath = fullfile(dataPath, fileName);
            
            % Double check file existence after string replacement
            if ~exist(dataFilePath, 'file')
                warning('File not found after cleaning name: %s', dataFilePath);
                continue;
            end
            
            %% Load data and calculate PCIst
            fprintf('Processing: Subj %s | %s | Epoch %d...\n', current_user, current_session, current_epoch);
            
            try
                % Load EEG variable into the workspace
                load(dataFilePath); 
                
                % Average across all trials (3rd dimension)
                tep_matrix = mean(EEG.data, 3);
                time_vector = EEG.times;
                
                % Execute PCIst algorithm
                [pci_value, dNST] = PCIst(tep_matrix, time_vector, params);
                
                % Store values into results table
                results.Subject{rowIdx}    = current_user;
                results.Session{rowIdx}    = current_session;
                results.Epoch(rowIdx)      = current_epoch;
                results.PCI_ST(rowIdx)     = pci_value;
                
                rowIdx = rowIdx + 1;
                
            catch ME
                warning('Error processing file %s: %s', fileName, ME.message);
            end
            
        end % End of Epoch loop
    end % End of Session loop
end % End of User loop

%% Display Final Summary Table
fprintf('\n=======================================================\n');
fprintf('                 PCI_ST PROCESSING COMPLETE            \n');
fprintf('=======================================================\n');
disp(results);

% Optional: Save results to a CSV file
writetable(results, fullfile(baseDataPath, 'PCIst_compiled_results.csv'));