clearvars;
cd('C:\Users\nasak\OneDrive\Documents\MATLAB\For Luis')
eeglab_path = 'C:\Users\nasak\OneDrive\Documents\MATLAB\eeglab2026.0.0';
addpath(eeglab_path);
eeglab;
addpath(genpath('C:\Users\nasak\OneDrive\Documents\MATLAB\For Luis\discover-eeg-master'))

ICLabel = [
    NaN NaN; %Brain
    0.8 1.0; %Muscle
    0.8 1.0; %Eye
    NaN NaN; %Heart
    NaN NaN; %Line Noise
    NaN NaN; %Channel Noise
    NaN NaN; %Other
];



%Location of the main study directory
DIR = 'C:\Users\nasak\OneDrive\Documents\MATLAB\For Luis\data'; %using source_data/Formatted instead of raw_data - RC
OutputBase = 'C:\Users\nasak\OneDrive\Documents\MATLAB\For Luis\data\outputs\resting_preprocessing';

%Location of the folder that contains this script and any associated processing files
Current_File_Path = 'C:\Users\nasak\OneDrive\Documents\MATLAB\For Luis';
ChanLocsDir='C:\Users\nasak\OneDrive\Documents\MATLAB\eeglab2026.0.0\plugins\dipfit\standard_BEM\elec\standard_1005.elc';

%Open EEGLAB and ERPLAB Toolboxes
cd(eeglab_path);
eeglab;

subs = dir(DIR);
subs = subs([subs.isdir]);
subs = subs(~ismember({subs.name}, {'.', '..', 'outputs'}));

sublist = 1:length(subs);

for xloop = 1:length(sublist)
    subject=sublist(xloop);

    subjectID = subs(subject).name;
    subjectfolder=fullfile(subs(subject).folder, subs(subject).name);

    dataFolder = fullfile(subjectfolder, 'data');
    if ~exist(dataFolder, 'dir')
        disp(['Skipping (no Data folder): ' subjectID]);
        continue;
    end

    ses = dir(dataFolder);
    ses = ses([ses.isdir]);
    ses = ses(~ismember({ses.name}, {'.', '..'}));
    seslist = 1:length(ses);

    for session=1:length(seslist)

        sessionName = ses(session).name;
        sessionFolder = fullfile(ses(session).folder, ses(session).name);
        
        sessionfolderpath = fullfile(sessionFolder, 'Resting EEG');
        if ~exist(sessionfolderpath, 'dir')
            continue;
        end


        vhdrFiles = dir(fullfile(sessionfolderpath, '*.vhdr'));

        vhdrFiles = vhdrFiles(~startsWith({vhdrFiles.name}, '._'));
        vhdrFiles = vhdrFiles(~startsWith({vhdrFiles.name}, '.'));


        if isempty(vhdrFiles)
            continue;
        end

        
        Outputfolderpath = fullfile(OutputBase, subjectID, sessionName);
        if ~exist(Outputfolderpath, 'dir')
            mkdir(Outputfolderpath)
        end

        for f= 1:length(vhdrFiles)

            %Define subject path based on study directory and subject ID of current subject
            % Subject_Path = [ses(j).folder filesep ses(j).name];
            vhdrName = vhdrFiles(f).name;
            vhdrStem = erase(vhdrName, '.vhdr');
            
            %skip redundanciesfolderpath
            postSetPath = fullfile(Outputfolderpath, ['PostICA_' subjectID '_' sessionName '_' vhdrStem '.set']);
            if exist(postSetPath, 'file')
                disp(['Skipping (already processed): ' subjectID ' / ' sessionName ' / ' vhdrName]);
                continue;
            end

            disp(['Processing: ' subjectID ' / ' sessionName ' / ' vhdrName]);

            %% STEP 1 % import data

            EEG = pop_loadbv(sessionfolderpath, vhdrName); %to use vhdr files - error with line 85 of biosig2eeglab (line 235 of pop_biosig)


            %% Step 3: look up channel locations

            % look up channel info

            EEG=pop_chanedit(EEG, 'lookup',ChanLocsDir); % need to be updated to include the actual function that tells the reference channels
            EEG=pop_select(EEG,'channel',1:min(70, EEG.nbchan));
            OGEEG=EEG;

            %% Step 4: Downsample

            EEG = pop_resample(EEG, 512);

            %% Step 5: Clean line noise
            EEG = pop_cleanline(EEG, 'linefreqs', 60, 'newversion', 0); % check if the newer version works with parallel processing parameters
            %% Step 6: Remove bad channels and Interpolate them

            EEG = pop_clean_rawdata(EEG,...
                'FlatlineCriterion',5,...
                'ChannelCriterion',0.8,...
                'LineNoiseCriterion',4,...
                'Highpass',[0.25 0.75],...
                'BurstCriterion',20,...
                'WindowCriterion',.25, ...
                'WindowCriterionTolerances', [-Inf,7], ...
                'BurstRejection','off',...
                'Distance','Euclidian');

            EEG = pop_interp(EEG,eeg_mergelocs(OGEEG.chanlocs),'spherical');
            %% Step 7: Rereference data
            % EEG = pop_reref( EEG, []);
            
            if EEG.nbchan >= 70
                EEG = pop_reref( EEG, [69, 70]);
            else
                EEG = pop_reref(EEG, []);
            end

            PREstudy = pop_saveset(EEG,'filename', ['PreICA_' subjectID '_' sessionName,'_' vhdrStem], 'filepath', Outputfolderpath);

            %% Step 8: Run ICA])])
            EEGtemp=EEG;
            EEGtemp = pop_runica(EEGtemp,'icatype','runica','concatcond','off');
            EEGtemp = pop_iclabel(EEGtemp,'default');
            EEGtemp = pop_icflag(EEGtemp, ICLabel); % flag artifactual components using IClabel

            classifications = EEGtemp.etc.ic_classification.ICLabel.classifications; % Keep classifications before component substraction
            EEGtemp = pop_subcomp(EEGtemp,[],0); % Subtract artifactual independent components
            EEGtemp.etc.ic_classification.ICLabel.orig_classifications = classifications;

            POSTstudy = pop_saveset(EEGtemp,'filename', ['PostICA_' subjectID '_' sessionName, '_' vhdrStem], 'filepath',Outputfolderpath);

        end
    end
end

