clear;
close all;
clc;

gsPath = '/usr/local/bin/gs';            % full gs executable
gsDir  = fileparts(gsPath);                 % -> '/opt/homebrew/bin'
currentPath = getenv('PATH');
setenv('PATH', [gsDir ':' currentPath]);    % prepend gs directory
setenv('GSCMD', gsPath);   

%% Setup Environment
% Define subject number and set the directory

usr_num = '005';

FileDir = '/Volumes/T7/OPTE';
% FileDir='/data/OPTE/Data';
newDir=[FileDir '/' usr_num '/Data' ];

% Define the session names 
sessions = {'Baseline TMSEEG', 'Dosing Session', '1-week follow-up', '2-week follow-up', '1-month follow-up'};


% =========================================================================
% LOOP THROUGH ALL SESSIONS
% =========================================================================
for session_num = 1:length(sessions)
    
    % Determine which epochs to cycle through based on the current session
    if strcmp(sessions{session_num}, 'Dosing Session')
        epochs_to_run = 1:4; % Cycle through epochs 1 to 4 for Dosing
    else
        epochs_to_run = 0;   % Only run epoch 0 for all other sessions
    end
    
    % LOOP THROUGH DETERMINED EPOCHS
    for epoch_num = epochs_to_run
        
        % Dynamic Protocol Mapping based on Epoch Number
        if epoch_num == 0 || epoch_num == 1
            protocol_num = 4;
        elseif epoch_num == 2
            protocol_num = 8;
        elseif epoch_num == 3
            protocol_num = 12;
        elseif epoch_num == 4
            protocol_num = 16;
        end
        
        % Reset and construct the file directory for the current session
        cd([newDir '/' sessions{session_num}])
        FileDir=[]; FileDir=[newDir '/' sessions{session_num} '/TMSEEG'];
        cd(FileDir)
        listtemp=dir('*vhdr'); % list all conditions within a directory
        listtemp = listtemp(~startsWith({listtemp.name}, '.')); % exclude files starting with "."
        prot_list={listtemp.name}';

        %% Load Data
        % Load the selected EEG data file as an EEGLAB structure
        [EEG, misc] = c_TMSEEG_prepareForPreprocessing(...
            'inputFilePath', [FileDir '/' prot_list{protocol_num}],...
            'pulseEvent', 'T  1', ... % Specify the TMS pulse event marker
            'epochTimespan', [-0.5 0.5]); % epoch time span in seconds

        %% Call the main preprocessing pipeline script
        EEG = c_TMSEEG_Preprocess_AARATEPPipeline(EEG,...
            'pulseEvent', misc.pulseEvent,...
            'epochTimespan', misc.epochTimespan,...
            'outputDir', ['/Volumes/T7/OPTE/' usr_num '/Preprocessed TMS/' sessions{session_num} '/Epoch' num2str(epoch_num)],...
            'outputFilePrefix', [usr_num '_' sessions{session_num} '_' prot_list{protocol_num}]);
            
    end % End of epoch loop
end % End of session loop
