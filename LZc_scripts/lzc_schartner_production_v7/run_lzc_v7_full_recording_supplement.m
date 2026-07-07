%% run_lzc_v7_full_recording_supplement.m
% Add full-recording Schartner-style LZ outputs to an existing V7 run.
%
% This script does NOT rerun the windowed V7 pipeline. It loads each
% PostICA .set file once and computes:
%   - full-recording channel-wise LZs/LZsN
%   - full-recording all-channel LZc/LZcN
%
% Outputs:
%   LZc_Schartner_ProductionV7_Results/CSV/
%       PUBLIC_full_recording_LZs_by_channel.csv
%       PUBLIC_full_recording_LZc_by_file.csv

clear; clc;
format long g

%% ---------------- PATH CONFIGURATION ----------------
% Set USE_HARDCODED_PATHS=true for reproducible batch runs.
% Leave it false to select input folders interactively when the script runs.
USE_HARDCODED_PATHS = false;

% Optional hardcoded paths. Edit these for your machine if needed.
HARD_CODED_EEGLAB_PATH = '';           % e.g., '/path/to/eeglab2026.0.0'
HARD_CODED_INPUT_ROOT = '';            % folder containing PostICA_*.set files
HARD_CODED_OUTPUT_ROOT = '';           % optional; leave blank to write inside input_root
HARD_CODED_V7_CODE_DIR = '';           % optional; defaults to this script folder
HARD_CODED_LEGACY_LZC_HELPER_DIR = ''; % optional helper folder if your setup keeps helpers elsewhere

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, '-begin');

eeglab_path = char(HARD_CODED_EEGLAB_PATH);
if isempty(char(HARD_CODED_V7_CODE_DIR))
    v7_code_dir = script_dir;
else
    v7_code_dir = char(HARD_CODED_V7_CODE_DIR);
end
legacy_lzc_scripts_dir = char(HARD_CODED_LEGACY_LZC_HELPER_DIR);

%% ---------------- SETTINGS ----------------

seed = 0;

% Match V7 defaults unless you intentionally want the 10-shuffle extension.
nBinaryShuffles = 1;
nPhaseSurrogates = 10;

overwrite_existing = false;

% Leave "" for all files.
target_participant = "";  % "005"
target_session = "";      % "baseline"
target_eyes = "";         % "EO"
target_epoch = "";        % "Epoch0"

%% ---------------- SETUP ----------------

add_optional_path(eeglab_path, 'EEGLAB');
add_optional_path(v7_code_dir, 'V7 code folder');
add_optional_path(legacy_lzc_scripts_dir, 'legacy LZc helper folder');

eeglab;

old_fig_visibility = get(groot, 'defaultFigureVisible');
set(groot, 'defaultFigureVisible', 'off');
cleanup_obj = onCleanup(@() set(groot, 'defaultFigureVisible', old_fig_visibility));

if USE_HARDCODED_PATHS
    input_root = char(HARD_CODED_INPUT_ROOT);
    if isempty(input_root)
        error('HARD_CODED_INPUT_ROOT is empty. Set it or use interactive folder selection.');
    end
else
    input_root = uigetdir(pwd, 'Select folder containing PostICA EEGLAB .set files');
    if isequal(input_root, 0)
        disp('Folder selection cancelled.');
        return;
    end
end

if ~exist(input_root, 'dir')
    error('Input folder not found: %s', input_root);
end

if ~isempty(char(HARD_CODED_OUTPUT_ROOT))
    output_root = char(HARD_CODED_OUTPUT_ROOT);
else
    output_root = fullfile(input_root, 'LZc_Schartner_ProductionV7_Results');
end
csv_root = fullfile(output_root, 'CSV');
per_file_root = fullfile(output_root, 'Per_File_FullRecording');

mkdir_if_needed(output_root);
mkdir_if_needed(csv_root);
mkdir_if_needed(per_file_root);

channel_csv = fullfile(csv_root, 'PUBLIC_full_recording_LZs_by_channel.csv');
file_csv = fullfile(csv_root, 'PUBLIC_full_recording_LZc_by_file.csv');

if overwrite_existing
    delete_if_exists(channel_csv);
    delete_if_exists(file_csv);
end

existing_channel = read_existing_table(channel_csv);
existing_file = read_existing_table(file_csv);

existing_file_keys = strings(0, 1);
if ~isempty(existing_file)
    existing_file_keys = make_file_keys(existing_file);
end

%% ---------------- FIND FILES ----------------

set_files = dir(fullfile(input_root, '**', 'PostICA_*.set'));

if isempty(set_files)
    error('No PostICA_*.set files found under: %s', input_root);
end

fprintf('\nFound %d PostICA .set files.\n', numel(set_files));

all_channel_rows = table();
all_file_rows = table();

%% ---------------- MAIN LOOP ----------------

for f = 1:numel(set_files)

    file = set_files(f).name;
    folder = set_files(f).folder;
    set_path = fullfile(folder, file);
    [~, base_name, ~] = fileparts(file);

    info = parse_lzc_filename(base_name);

    if info.skip
        fprintf('Skipping unparsable file: %s\n', file);
        continue;
    end

    info.participant = normalize_participant_id(info.participant);

    if ~matches_targets(info, target_participant, target_session, target_eyes, target_epoch)
        continue;
    end

    file_key = make_one_file_key(info.participant, info.session, info.eyes, info.epoch, string(base_name));

    if ~overwrite_existing && any(existing_file_keys == file_key)
        fprintf('Skipping existing full-recording output: %s\n', file_key);
        continue;
    end

    fprintf('\n====================================================\n');
    fprintf('Full-recording supplement %d/%d:\n%s\n', f, numel(set_files), set_path);

    EEG = pop_loadset('filename', file, 'filepath', folder);
    X = double(EEG.data);
    [nChannels, nSamples] = size(X);
    srate = double(EEG.srate);
    chan_labels = get_channel_labels(EEG, nChannels);

    participant_dir = fullfile(per_file_root, char(info.participant));
    mkdir_if_needed(participant_dir);

    file_tag = sprintf('%s_%s_%s_%s_%s', ...
        char(info.participant), char(info.session), char(info.eyes), char(info.epoch), base_name);

    fprintf('Computing full-recording channel-wise LZs/LZsN for %d channels.\n', nChannels);

    channel_rows = table();

    for ch = 1:nChannels

        fprintf('  Channel %d/%d: %s\n', ch, nChannels, chan_labels(ch));

        ch_seed = seed + ch;
        lzs = v5_lzs_segment(X(ch, :), ch_seed, nBinaryShuffles, nPhaseSurrogates);

        row = table( ...
            string(info.participant), string(info.session), string(info.eyes), string(info.epoch), ...
            string(base_name), string(set_path), ...
            chan_labels(ch), double(ch), ...
            double(lzs.LZs), double(lzs.LZsN), double(lzs.RawLZs), ...
            double(lzs.BinaryShuffleMeanRawLZs), double(lzs.BinaryShuffleSDRawLZs), ...
            double(lzs.PhaseRawLZsMean), double(lzs.PhaseRawLZsSD), ...
            double(lzs.Threshold), double(lzs.PropOnes), double(lzs.NTransitions), ...
            double(lzs.StringLength), double(nSamples), double(srate), ...
            double(nBinaryShuffles), double(nPhaseSurrogates), ...
            "FullRecordingV7Supplement", ...
            'VariableNames', {'Participant','Session','Eyes','Epoch','File','Path', ...
            'Channel','ChannelIndex','LZs','LZsN','RawLZs', ...
            'BinaryShuffleMeanRawLZs','BinaryShuffleSDRawLZs', ...
            'PhaseRawLZsMean','PhaseRawLZsSD','Threshold','PropOnes','NTransitions', ...
            'StringLength','N_Samples','SamplingRate','N_BinaryShuffles','N_PhaseSurrogates', ...
            'PipelineVersion'} ...
        );

        channel_rows = [channel_rows; row];

    end

    fprintf('Computing full-recording all-channel LZc/LZcN.\n');

    lzc = v5_lzc_all_channels_segment(X, seed, nBinaryShuffles, nPhaseSurrogates);

    file_row = table( ...
        string(info.participant), string(info.session), string(info.eyes), string(info.epoch), ...
        string(base_name), string(set_path), ...
        double(lzc.LZc), double(lzc.LZcN), double(lzc.RawLZc), ...
        double(lzc.BinaryShuffleMeanRawLZc), double(lzc.BinaryShuffleSDRawLZc), ...
        double(lzc.PhaseRawLZcMean), double(lzc.PhaseRawLZcSD), ...
        double(lzc.PropOnes), double(lzc.NTransitions), double(lzc.StringLength), ...
        double(nChannels), double(nSamples), double(srate), ...
        double(nBinaryShuffles), double(nPhaseSurrogates), ...
        "FullRecordingV7Supplement", ...
        'VariableNames', {'Participant','Session','Eyes','Epoch','File','Path', ...
        'LZc','LZcN','RawLZc','BinaryShuffleMeanRawLZc','BinaryShuffleSDRawLZc', ...
        'PhaseRawLZcMean','PhaseRawLZcSD','PropOnes','NTransitions','StringLength', ...
        'N_Channels','N_Samples','SamplingRate','N_BinaryShuffles','N_PhaseSurrogates', ...
        'PipelineVersion'} ...
    );

    per_channel_path = fullfile(participant_dir, [file_tag '_full_recording_LZs_by_channel.csv']);
    per_file_path = fullfile(participant_dir, [file_tag '_full_recording_LZc_by_file.csv']);

    writetable(channel_rows, per_channel_path);
    writetable(file_row, per_file_path);

    all_channel_rows = [all_channel_rows; channel_rows];
    all_file_rows = [all_file_rows; file_row];

    fprintf('Saved per-file full-recording supplement:\n%s\n%s\n', per_channel_path, per_file_path);

end

%% ---------------- SAVE PUBLIC CSVs ----------------

if ~isempty(all_channel_rows)
    if ~isempty(existing_channel)
        all_channel_rows = [existing_channel; all_channel_rows];
    end
    all_channel_rows = sortrows(all_channel_rows, {'Participant','Session','Eyes','Epoch','File','ChannelIndex'});
    writetable(all_channel_rows, channel_csv);
    fprintf('\nSaved full-recording channel table:\n%s\n', channel_csv);
else
    fprintf('\nNo new full-recording channel rows were created.\n');
end

if ~isempty(all_file_rows)
    if ~isempty(existing_file)
        all_file_rows = [existing_file; all_file_rows];
    end
    all_file_rows = sortrows(all_file_rows, {'Participant','Session','Eyes','Epoch','File'});
    writetable(all_file_rows, file_csv);
    fprintf('Saved full-recording file table:\n%s\n', file_csv);
else
    fprintf('No new full-recording file rows were created.\n');
end

fprintf('\nFull-recording supplement complete.\n');

%% ---------------- LOCAL FUNCTIONS ----------------

function add_optional_path(p, label)
p = char(p);
if isempty(p)
    fprintf('%s path not set; continuing if required functions are already on the MATLAB path.\n', label);
    return;
end
if ~exist(p, 'dir')
    error('%s not found: %s', label, p);
end
addpath(p, '-begin');
end

function mkdir_if_needed(p)
if ~exist(p, 'dir')
    mkdir(p);
end
end

function delete_if_exists(p)
if exist(p, 'file')
    delete(p);
end
end

function T = read_existing_table(p)
if exist(p, 'file')
    T = readtable(p, 'TextType', 'string');
else
    T = table();
end
end

function tf = matches_targets(info, target_participant, target_session, target_eyes, target_epoch)
tf = true;
if strlength(string(target_participant)) > 0
    tf = tf && string(info.participant) == normalize_participant_id(target_participant);
end
if strlength(string(target_session)) > 0
    tf = tf && string(info.session) == string(target_session);
end
if strlength(string(target_eyes)) > 0
    tf = tf && string(info.eyes) == string(target_eyes);
end
if strlength(string(target_epoch)) > 0
    tf = tf && string(info.epoch) == string(target_epoch);
end
end

function keys = make_file_keys(T)
keys = strings(height(T), 1);
for i = 1:height(T)
    keys(i) = make_one_file_key(T.Participant(i), T.Session(i), T.Eyes(i), T.Epoch(i), T.File(i));
end
end

function key = make_one_file_key(participant, session_name, eyes, epoch, file)
key = normalize_participant_id(participant) + "|" + string(session_name) + "|" + ...
    string(eyes) + "|" + string(epoch) + "|" + string(file);
end
