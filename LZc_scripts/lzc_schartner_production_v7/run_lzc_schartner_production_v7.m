%% run_lzc_schartner_production_v7.m
% Production V7 Schartner-style LZs/LZc pipeline.
%
% Purpose of Production V7:
%   Fix the summary/bookkeeping issue once and for all.
%
% Design rule:
%   - Window-level rows are the source of truth.
%   - File summaries are recomputed ONLY from window-level numeric columns.
%   - Summary values are never trusted if they cannot be reconstructed.
%   - Validation checks stop the script immediately if any ratio is wrong.
%
% Method:
%   LZs  = RawLZs / BinaryShuffleMeanRawLZs
%   LZsN = RawLZs / PhaseRawLZsMean
%
%   LZc  = RawLZc / BinaryShuffleMeanRawLZc
%   LZcN = RawLZc / PhaseRawLZcMean
%
% Schartner core:
%   mean-center + detrend -> Hilbert envelope -> mean threshold -> binary string
%   -> Schartner cpr dictionary count.
%
% Main V7 branch uses all channels for LZc. Use ../paper_method_random10 for Schartner random 10-channel picks.

clear; clc;
format long g

%% ---------------- PATH CONFIGURATION ----------------
% Set USE_HARDCODED_PATHS=true for reproducible batch runs.
% Leave it false to select input folders interactively when the script runs.
USE_HARDCODED_PATHS = false;

% Optional hardcoded paths. Edit these for your machine if USE_HARDCODED_PATHS=true.
HARD_CODED_EEGLAB_PATH = '';      % e.g., '/path/to/eeglab2026.0.0'
HARD_CODED_INPUT_ROOT = '';       % folder containing PostICA_*.set files
HARD_CODED_OUTPUT_ROOT = '';      % optional; leave blank to write inside input_root

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, '-begin');

eeglab_path = char(HARD_CODED_EEGLAB_PATH);
if ~isempty(eeglab_path)
    if ~exist(eeglab_path, 'dir')
        error('EEGLAB path not found: %s', eeglab_path);
    end
    addpath(eeglab_path, '-begin');
elseif exist('eeglab', 'file') ~= 2
    warning(['EEGLAB path is blank and eeglab is not currently on the MATLAB path. ' ...
        'Set HARD_CODED_EEGLAB_PATH if eeglab does not start.']);
end

eeglab;

%% ---------------- SETTINGS ----------------

seed = 0;

window_lengths_sec = [2 10];
window_overlap_sec = 0;

% Exact original Schartner code uses one binary shuffle.
% Use 1 for replication; use 10 only as an intentional stability extension.
nBinaryShuffles = 1;

% Use 1 for quick testing, 10 for main run.
nPhaseSurrogates = 10;

overwrite_existing = false;

save_channel_long_csv = true;

% Long-format accessibility outputs:
%   true  = save all channel-window rows and all LZc-window rows.
%   false = only save summaries.
save_long_csvs_for_sharing = true;

% Optional debug outputs. These can become large.
% Main run recommendation:
%   save_debug_mat = false;
% For one-file debugging:
%   save_debug_mat = true;
save_debug_mat = false;

% When save_debug_mat=true:
%   "compact" saves metrics/thresholds/binary-string snippets only.
%   "full" also saves B, TH, M for selected windows/channels and can be large.
debug_mode = "compact";

% Only save full debug arrays for the first N windows/channels per file-window size.
% This prevents accidental enormous .mat files.
debug_max_windows = 2;
debug_max_channels = 5;

% Processing provenance metadata written into summary outputs.
pipeline_version = "ProductionV7";
date_processed = string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss"));

% Quick-test filters. Leave "" for all.
target_participant = "";  % "005"
target_session = "";      % "baseline"
target_eyes = "";         % "EC"
target_epoch = "";        % "Epoch0"

% Figure suppression
old_fig_visibility = get(groot, 'defaultFigureVisible');
set(groot, 'defaultFigureVisible', 'off');
cleanup_obj = onCleanup(@() set(groot, 'defaultFigureVisible', old_fig_visibility));

%% ---------------- INPUT / OUTPUT ----------------

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
per_file_root = fullfile(output_root, 'Per_File');
debug_root = fullfile(output_root, 'Debug_MAT');

v5_mkdir(output_root);
v5_mkdir(csv_root);
v5_mkdir(per_file_root);
v5_mkdir(debug_root);

paths.parser = fullfile(csv_root, 'parser_check.csv');
paths.file_summary = fullfile(csv_root, 'ALL_file_summary.csv');
paths.window_summary = fullfile(csv_root, 'ALL_window_summary.csv');
paths.channel_long = fullfile(csv_root, 'ALL_channel_LZs_long.csv');
paths.channel_summary = fullfile(csv_root, 'ALL_channel_LZs_summary.csv');
paths.lzc_long = fullfile(csv_root, 'ALL_all_channel_LZc_long.csv');
paths.window_comparison = fullfile(csv_root, 'ALL_window_size_comparison.csv');

% Extra public-facing long-format files with intentionally plain names.
paths.public_window_channel_lzs = fullfile(csv_root, 'PUBLIC_long_LZs_by_channel_window.csv');
paths.public_window_lzc = fullfile(csv_root, 'PUBLIC_long_LZc_by_window.csv');
paths.public_file_summary = fullfile(csv_root, 'PUBLIC_file_summary.csv');


if overwrite_existing
    v5_delete(paths.parser);
    v5_delete(paths.file_summary);
    v5_delete(paths.window_summary);
    v5_delete(paths.channel_long);
    v5_delete(paths.channel_summary);
    v5_delete(paths.lzc_long);
    v5_delete(paths.window_comparison);
    v5_delete(paths.public_window_channel_lzs);
    v5_delete(paths.public_window_lzc);
    v5_delete(paths.public_file_summary);
end

%% ---------------- FIND FILES ----------------

set_files = dir(fullfile(input_root, '**', 'PostICA_*.set'));

if isempty(set_files)
    error('No PostICA_*.set files found under: %s', input_root);
end

fprintf('\nFound %d PostICA .set files.\n', numel(set_files));

%% ---------------- PARSE / FILTER ----------------

parser_table = table();
keep = false(numel(set_files), 1);

for f = 1:numel(set_files)

    file = set_files(f).name;
    folder = set_files(f).folder;
    set_path = fullfile(folder, file);
    [~, base_name, ~] = fileparts(file);

    info = parse_lzc_filename(base_name);

    if info.skip
        continue;
    end

    info.participant = normalize_participant_id(info.participant);

    row = table( ...
        string(info.participant), string(info.session), string(info.eyes), string(info.epoch), ...
        string(base_name), string(set_path), ...
        'VariableNames', {'Participant','Session','Eyes','Epoch','File','Path'} ...
    );

    parser_table = [parser_table; row];

    matches = true;

    if strlength(string(target_participant)) > 0
        matches = matches && string(info.participant) == normalize_participant_id(target_participant);
    end
    if strlength(string(target_session)) > 0
        matches = matches && string(info.session) == string(target_session);
    end
    if strlength(string(target_eyes)) > 0
        matches = matches && string(info.eyes) == string(target_eyes);
    end
    if strlength(string(target_epoch)) > 0
        matches = matches && string(info.epoch) == string(target_epoch);
    end

    keep(f) = matches;

end

parser_table = normalize_lzc_table_ids(parser_table);
writetable(parser_table, paths.parser);

fprintf('\nParser summary by participant:\n');
if ~isempty(parser_table)
    disp(groupsummary(parser_table, 'Participant'));
end

fprintf('\nParser summary by session/eyes/epoch:\n');
if ~isempty(parser_table)
    disp(groupsummary(parser_table, {'Session','Eyes','Epoch'}));
end

set_files = set_files(keep);

if isempty(set_files)
    error('No files matched target filters.');
end

fprintf('\nRunning Production V7 Schartner pipeline on %d files.\n', numel(set_files));

%% ---------------- MAIN LOOP ----------------

for f = 1:numel(set_files)

    file = set_files(f).name;
    folder = set_files(f).folder;
    set_path = fullfile(folder, file);
    [~, base_name, ~] = fileparts(file);

    info = parse_lzc_filename(base_name);
    if info.skip
        continue;
    end
    info.participant = normalize_participant_id(info.participant);

    file_tag = sprintf('%s_%s_%s_%s_%s', ...
        char(info.participant), char(info.session), char(info.eyes), char(info.epoch), base_name);

    participant_dir = fullfile(per_file_root, char(info.participant));
    v5_mkdir(participant_dir);

    done_marker = fullfile(participant_dir, [file_tag '_DONE.txt']);

    if ~overwrite_existing && exist(done_marker, 'file')
        fprintf('\nSkipping completed file: %s\n', file_tag);
        continue;
    end

    fprintf('\n====================================================\n');
    fprintf('Processing %d/%d:\n%s\n', f, numel(set_files), set_path);

    EEG = pop_loadset('filename', file, 'filepath', folder);
    X = double(EEG.data);
    [nChannels, nSamples] = size(X);
    srate = double(EEG.srate);
    chan_labels = get_channel_labels(EEG, nChannels);

    per_file_window_summaries = table();

    for w = 1:numel(window_lengths_sec)

        window_sec = double(window_lengths_sec(w));
        window_samples = round(window_sec * srate);
        overlap_samples = round(window_overlap_sec * srate);

        windows = v5_make_windows(nSamples, window_samples, overlap_samples, srate);

        if isempty(windows)
            warning('No valid windows for %s at %.3f s.', file_tag, window_sec);
            continue;
        end

        fprintf('\nWindow %.3f sec: %d windows\n', window_sec, height(windows));

        window_tag = sprintf('%s_window%gs', file_tag, window_sec);

        per_window_summary_path = fullfile(participant_dir, [window_tag '_window_summary.csv']);
        per_channel_long_path = fullfile(participant_dir, [window_tag '_channel_LZs_long.csv']);
        per_channel_summary_path = fullfile(participant_dir, [window_tag '_channel_LZs_summary.csv']);
        per_lzc_long_path = fullfile(participant_dir, [window_tag '_all_channel_LZc_long.csv']);
        per_file_summary_path = fullfile(participant_dir, [window_tag '_file_summary.csv']);

        if ~overwrite_existing && exist(per_file_summary_path, 'file')
            fprintf('Window result exists. Loading: %s\n', per_file_summary_path);
            T = readtable(per_file_summary_path, 'TextType', 'string');
            per_file_window_summaries = [per_file_window_summaries; T];
            continue;
        end

        channel_long = table();
        lzc_long = table();
        debug_records = struct();
        debug_records.file_tag = file_tag;
        debug_records.window_sec = window_sec;
        debug_records.nBinaryShuffles = nBinaryShuffles;
        debug_records.nPhaseSurrogates = nPhaseSurrogates;
        debug_records.debug_mode = debug_mode;
        debug_records.lzc_windows = struct([]);
        debug_records.lzs_channels = struct([]);

        %% ---------------- WINDOW LOOP ----------------

        for seg = 1:height(windows)

            idx = windows.StartSample(seg):windows.EndSample(seg);

            fprintf('  Window %d/%d | %.3f-%.3f sec\n', ...
                seg, height(windows), windows.StartSec(seg), windows.EndSec(seg));

            %% LZs by channel

            for ch = 1:nChannels

                lzs = v5_lzs_segment( ...
                    X(ch, idx), seed + seg*1000 + ch, nBinaryShuffles, nPhaseSurrogates);

                v5_validate_identity("LZs", lzs.LZs, lzs.RawLZs, lzs.BinaryShuffleMeanRawLZs, ...
                    sprintf('%s window %d channel %s', file_tag, seg, chan_labels(ch)));

                v5_validate_identity("LZsN", lzs.LZsN, lzs.RawLZs, lzs.PhaseRawLZsMean, ...
                    sprintf('%s window %d channel %s', file_tag, seg, chan_labels(ch)));

                if save_debug_mat && seg <= debug_max_windows && ch <= debug_max_channels
                    drec = v5_make_debug_lzs_record(lzs, seg, windows.StartSec(seg), windows.EndSec(seg), ch, string(chan_labels(ch)), debug_mode);
                    debug_records.lzs_channels(end+1) = drec; %#ok<SAGROW>
                end

                row = table( ...
                    string(info.participant), string(info.session), string(info.eyes), string(info.epoch), string(base_name), ...
                    double(window_sec), double(seg), double(windows.StartSec(seg)), double(windows.EndSec(seg)), ...
                    string(chan_labels(ch)), double(ch), ...
                    double(lzs.LZs), double(lzs.LZsN), ...
                    double(lzs.RawLZs), double(lzs.BinaryShuffleMeanRawLZs), double(lzs.BinaryShuffleSDRawLZs), ...
                    double(lzs.PhaseRawLZsMean), double(lzs.PhaseRawLZsSD), ...
                    double(lzs.Threshold), double(lzs.PropOnes), double(lzs.NTransitions), double(lzs.StringLength), ...
                    'VariableNames', {'Participant','Session','Eyes','Epoch','File', ...
                    'WindowSec','Window','WindowStartSec','WindowEndSec', ...
                    'Channel','ChannelIndex', ...
                    'LZs','LZsN', ...
                    'RawLZs','BinaryShuffleMeanRawLZs','BinaryShuffleSDRawLZs', ...
                    'PhaseRawLZsMean','PhaseRawLZsSD', ...
                    'Threshold','PropOnes','NTransitions','StringLength'} ...
                );

                channel_long = [channel_long; row];

            end

            %% LZc all channels

            lzc = v5_lzc_all_channels_segment( ...
                X(:, idx), seed + 1000000 + seg, nBinaryShuffles, nPhaseSurrogates);

            v5_validate_identity("LZc", lzc.LZc, lzc.RawLZc, lzc.BinaryShuffleMeanRawLZc, ...
                sprintf('%s window %d', file_tag, seg));

            v5_validate_identity("LZcN", lzc.LZcN, lzc.RawLZc, lzc.PhaseRawLZcMean, ...
                sprintf('%s window %d', file_tag, seg));

            if save_debug_mat && seg <= debug_max_windows
                drec = v5_make_debug_lzc_record(lzc, seg, windows.StartSec(seg), windows.EndSec(seg), debug_mode);
                debug_records.lzc_windows(end+1) = drec; %#ok<SAGROW>
            end

            lzc_row = table( ...
                string(info.participant), string(info.session), string(info.eyes), string(info.epoch), string(base_name), ...
                double(window_sec), double(seg), double(windows.StartSec(seg)), double(windows.EndSec(seg)), ...
                double(lzc.LZc), double(lzc.LZcN), ...
                double(lzc.RawLZc), double(lzc.BinaryShuffleMeanRawLZc), double(lzc.BinaryShuffleSDRawLZc), ...
                double(lzc.PhaseRawLZcMean), double(lzc.PhaseRawLZcSD), ...
                double(lzc.PropOnes), double(lzc.NTransitions), double(lzc.StringLength), double(nChannels), ...
                'VariableNames', {'Participant','Session','Eyes','Epoch','File', ...
                'WindowSec','Window','WindowStartSec','WindowEndSec', ...
                'LZc','LZcN', ...
                'RawLZc','BinaryShuffleMeanRawLZc','BinaryShuffleSDRawLZc', ...
                'PhaseRawLZcMean','PhaseRawLZcSD', ...
                'PropOnes','NTransitions','StringLength','N_Channels'} ...
            );

            lzc_long = [lzc_long; lzc_row];

        end

        %% ---------------- BUILD SUMMARIES FROM LONG TABLES ONLY ----------------

        window_summary = v5_build_window_summary_from_long(channel_long, lzc_long);
        channel_summary = v5_build_channel_summary_from_long(channel_long);
        file_summary = v5_build_file_summary_from_window_summary(window_summary, nSamples, srate, nBinaryShuffles, nPhaseSurrogates);

        file_summary.ReferenceStyle = "SchartnerProductionV7";
        file_summary.PipelineVersion = pipeline_version;
        file_summary.DateProcessed = date_processed;

        %% ---------------- VALIDATE SUMMARIES ----------------

        v5_validate_window_summary_against_long(window_summary, channel_long, lzc_long);
        v5_validate_file_summary_against_window_summary(file_summary, window_summary);

        %% ---------------- NORMALIZE IDs AND WRITE ----------------

        channel_long = normalize_lzc_table_ids(channel_long);
        channel_summary = normalize_lzc_table_ids(channel_summary);
        lzc_long = normalize_lzc_table_ids(lzc_long);
        window_summary = normalize_lzc_table_ids(window_summary);
        file_summary = normalize_lzc_table_ids(file_summary);

        writetable(window_summary, per_window_summary_path);
        writetable(file_summary, per_file_summary_path);
        writetable(channel_summary, per_channel_summary_path);
        writetable(lzc_long, per_lzc_long_path);

        if save_channel_long_csv
            writetable(channel_long, per_channel_long_path);
        end

        append_csv_v5(window_summary, paths.window_summary);
        append_csv_v5(file_summary, paths.file_summary);
        append_csv_v5(file_summary, paths.window_comparison);
        append_csv_v5(channel_summary, paths.channel_summary);
        append_csv_v5(lzc_long, paths.lzc_long);

        if save_channel_long_csv
            append_csv_v5(channel_long, paths.channel_long);
        end

        if save_long_csvs_for_sharing
            append_csv_v5(v5_public_channel_lzs_table(channel_long), paths.public_window_channel_lzs);
            append_csv_v5(v5_public_lzc_table(lzc_long), paths.public_window_lzc);
            append_csv_v5(v5_public_file_summary_table(file_summary), paths.public_file_summary);
        end

        if save_debug_mat
            participant_debug_dir = fullfile(debug_root, char(info.participant));
            v5_mkdir(participant_debug_dir);
            debug_mat_path = fullfile(participant_debug_dir, [window_tag '_debug_records.mat']);
            save(debug_mat_path, 'debug_records', '-v7.3');
        end

        per_file_window_summaries = [per_file_window_summaries; file_summary];

    end

    per_file_window_summaries = normalize_lzc_table_ids(per_file_window_summaries);
    writetable(per_file_window_summaries, fullfile(participant_dir, [file_tag '_file_summary_all_window_sizes.csv']));

    fid = fopen(done_marker, 'w');
    fprintf(fid, 'Completed %s at %s\n', file_tag, datestr(now));
    fclose(fid);

end

fprintf('\nProduction V7 complete. CSV outputs saved here:\n%s\n', csv_root);

%% ---------------- LOCAL HELPERS ----------------

function v5_mkdir(path_in)
if ~exist(path_in, 'dir')
    mkdir(path_in);
end
end

function v5_delete(path_in)
if exist(path_in, 'file')
    delete(path_in);
end
end
