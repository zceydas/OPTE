%% lzc_paper_method_pipeline.m
% Paper-method LZc/LZcN pipeline using your current LZc helper functions.
%
% This follows the Schartner-style method more closely:
%   1. Split each recording into non-overlapping 2-second segments.
%   2. For each segment, compute ordinary binary-shuffle-normalized LZc.
%   3. Phase-randomize each channel's time series within that segment.
%   4. Recompute LZc on phase-randomized surrogate segment data.
%   5. Compute LZcN = real LZc / mean(phase-surrogate LZc).
%   6. Average segment-level values to get file/condition summaries.
%
% This script produces:
%   - LZs by channel/segment: single-channel temporal diversity
%   - LZc by random channel-pick/segment: multi-channel spatiotemporal diversity
%   - File-level summaries
%   - Combined summary CSVs across all selected files
%
% Required helper functions on MATLAB path:
%   parse_lzc_filename.m
%   normalize_participant_id.m
%   normalize_lzc_table_ids.m
%   get_channel_labels.m
%   LZc_baseline_multishuffle.m
%   phase_shuffle_signal.m
%   pre_lzc.m
%   str_col_lzc.m
%   cpr_lzc.m

clear; clc;
format long g

%% ---------------- PATH CONFIGURATION ----------------
% Set USE_HARDCODED_PATHS=true for reproducible batch runs.
% Leave it false to select input folders interactively when the script runs.
USE_HARDCODED_PATHS = false;

% Optional hardcoded paths. Edit these for your machine if USE_HARDCODED_PATHS=true.
HARD_CODED_EEGLAB_PATH = '';    % e.g., '/path/to/eeglab2026.0.0'
HARD_CODED_INPUT_ROOT = '';     % folder containing PostICA_*.set files
HARD_CODED_OUTPUT_ROOT = '';    % optional; leave blank to write inside input_root

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

%% ---------------- USER SETTINGS ----------------

seed = 0;

% Paper used 2-second non-overlapping segments.
segment_sec = 2;
segment_overlap_sec = 0;

% Keep your current ordinary LZc normalization stable.
nBinaryShuffles = 10;

% Phase-randomized surrogate repetitions.
% Set to 1 for quick debugging, 10 for the main run.
nPhaseSurrogates = 10;

% Paper-style multichannel LZc:
% "30 random picks of 10 channels" per segment.
do_LZc_random_channel_picks = true;
nChannelPicks = 30;
nChannelsPerPick = 10;

% Paper-style single-channel LZs:
% compute per channel, then average across channels.
do_LZs_by_channel = true;

% Output detail.
save_segment_tables = true;      % large but useful for debugging/statistics
save_pick_tables = true;         % larger; stores all 30 random picks per segment
save_channel_tables = true;      % stores all channel x segment LZs values

% Skip already-computed file outputs.
overwrite_existing = false;

% Optional quick test restriction.
% Set to "" to run all available participants.
target_participant = "";         % e.g. "005"
target_session = "";             % e.g. "baseline" or "dosing"
target_eyes = "";                % e.g. "EC" or "EO"
target_epoch = "";               % e.g. "Epoch0", "Epoch1"

%% ---------------- QUIET FIGURE MODE ----------------

old_default_figure_visible = get(groot, 'defaultFigureVisible');
set(groot, 'defaultFigureVisible', 'off');
cleanup_obj = onCleanup(@() set(groot, 'defaultFigureVisible', old_default_figure_visible));

%% ---------------- SELECT INPUT FOLDER ----------------

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
    output_root = fullfile(input_root, 'LZc_PaperMethod_Results');
end

if ~exist(output_root, 'dir')
    mkdir(output_root);
end

%% ---------------- FIND FILES ----------------

set_files = dir(fullfile(input_root, '**', 'PostICA_*.set'));

if isempty(set_files)
    error('No PostICA_*.set files found under: %s', input_root);
end

fprintf('\nFound %d PostICA .set files.\n', numel(set_files));

%% ---------------- PARSER CHECK AND FILTER ----------------

parser_check = table();

keep_file = false(numel(set_files), 1);

for f = 1:numel(set_files)

    [~, base_name, ~] = fileparts(set_files(f).name);
    info = parse_lzc_filename(base_name);

    if info.skip
        continue;
    end

    info.participant = normalize_participant_id(info.participant);

    row = table( ...
        string(info.participant), ...
        string(info.session), ...
        string(info.eyes), ...
        string(info.epoch), ...
        string(base_name), ...
        string(fullfile(set_files(f).folder, set_files(f).name)), ...
        'VariableNames', {'Participant','Session','Eyes','Epoch','File','Path'} ...
    );

    parser_check = [parser_check; row];

    matches_target = true;

    if strlength(string(target_participant)) > 0
        matches_target = matches_target && string(info.participant) == normalize_participant_id(target_participant);
    end

    if strlength(string(target_session)) > 0
        matches_target = matches_target && string(info.session) == string(target_session);
    end

    if strlength(string(target_eyes)) > 0
        matches_target = matches_target && string(info.eyes) == string(target_eyes);
    end

    if strlength(string(target_epoch)) > 0
        matches_target = matches_target && string(info.epoch) == string(target_epoch);
    end

    keep_file(f) = matches_target;

end

parser_check = normalize_lzc_table_ids(parser_check);

writetable(parser_check, fullfile(output_root, 'paper_method_parser_check.csv'));

fprintf('\nParser summary by participant:\n');
if ~isempty(parser_check)
    disp(groupsummary(parser_check, 'Participant'));
end

fprintf('\nParser summary by session/eyes/epoch:\n');
if ~isempty(parser_check)
    disp(groupsummary(parser_check, {'Session','Eyes','Epoch'}));
end

set_files = set_files(keep_file);

if isempty(set_files)
    error('No files matched the target filters.');
end

fprintf('\nRunning paper-method pipeline on %d files after filtering.\n', numel(set_files));

%% ---------------- MAIN LOOP ----------------

all_file_summary = table();
all_LZs_channel_summary = table();
all_LZc_segment_summary = table();

for f = 1:numel(set_files)

    file = set_files(f).name;
    path = set_files(f).folder;
    set_path = fullfile(path, file);
    [~, base_name, ~] = fileparts(file);

    info = parse_lzc_filename(base_name);

    if info.skip
        fprintf('Skipping unparseable file: %s\n', set_path);
        continue;
    end

    info.participant = normalize_participant_id(info.participant);

    participant_dir = fullfile(output_root, char(info.participant));

    if ~exist(participant_dir, 'dir')
        mkdir(participant_dir);
    end

    file_tag = sprintf('%s_%s_%s_%s_%s', ...
        char(info.participant), char(info.session), char(info.eyes), char(info.epoch), base_name);

    summary_csv = fullfile(participant_dir, [file_tag '_paper_method_file_summary.csv']);
    LZs_channel_csv = fullfile(participant_dir, [file_tag '_paper_method_LZs_channel_summary.csv']);
    LZs_segment_csv = fullfile(participant_dir, [file_tag '_paper_method_LZs_channel_segments.csv']);
    LZc_segment_csv = fullfile(participant_dir, [file_tag '_paper_method_LZc_segment_summary.csv']);
    LZc_pick_csv = fullfile(participant_dir, [file_tag '_paper_method_LZc_random_picks.csv']);
    mat_path = fullfile(participant_dir, [file_tag '_paper_method_outputs.mat']);

    fprintf('\n====================================================\n');
    fprintf('Processing %d/%d:\n%s\n', f, numel(set_files), set_path);

    if ~overwrite_existing && exist(summary_csv, 'file') && exist(mat_path, 'file')

        fprintf('Paper-method outputs already exist. Loading summaries:\n%s\n', summary_csv);

        file_summary = readtable(summary_csv, 'TextType', 'string');
        file_summary = normalize_lzc_table_ids(file_summary);

        if exist(LZs_channel_csv, 'file')
            LZs_channel_summary = readtable(LZs_channel_csv, 'TextType', 'string');
            LZs_channel_summary = normalize_lzc_table_ids(LZs_channel_summary);
        else
            LZs_channel_summary = table();
        end

        if exist(LZc_segment_csv, 'file')
            LZc_segment_summary = readtable(LZc_segment_csv, 'TextType', 'string');
            LZc_segment_summary = normalize_lzc_table_ids(LZc_segment_summary);
        else
            LZc_segment_summary = table();
        end

        all_file_summary = [all_file_summary; file_summary];
        all_LZs_channel_summary = [all_LZs_channel_summary; LZs_channel_summary];
        all_LZc_segment_summary = [all_LZc_segment_summary; LZc_segment_summary];

        continue;

    end

    %% ---------------- LOAD EEG ----------------

    EEG = pop_loadset('filename', file, 'filepath', path);
    X = double(EEG.data);
    [nChannels, nSamples] = size(X);
    srate = double(EEG.srate);
    chan_labels = get_channel_labels(EEG, nChannels);

    segment_samples = round(segment_sec * srate);
    step_samples = round((segment_sec - segment_overlap_sec) * srate);

    if step_samples <= 0
        error('segment_overlap_sec must be less than segment_sec.');
    end

    seg_starts = 1:step_samples:(nSamples - segment_samples + 1);
    nSegments = numel(seg_starts);

    if nSegments < 1
        warning('File too short for one segment. Skipping: %s', set_path);
        continue;
    end

    fprintf('Channels: %d | Samples: %d | Segments: %d of %.3f sec\n', ...
        nChannels, nSamples, nSegments, segment_sec);

    %% ---------------- LZs: SINGLE-CHANNEL SEGMENT VALUES ----------------

    LZs_segment_table = table();
    LZs_channel_summary = table();

    if do_LZs_by_channel

        fprintf('\nComputing LZs by channel and 2-sec segment...\n');

        LZs_segment_table = table();

        for seg = 1:nSegments

            idx = seg_starts(seg):(seg_starts(seg) + segment_samples - 1);
            segment_start_sec = (seg_starts(seg) - 1) / srate;
            segment_end_sec = segment_start_sec + segment_sec;

            fprintf('  LZs segment %d/%d\n', seg, nSegments);

            for ch = 1:nChannels

                x_seg = X(ch, idx);

                [real_LZs, real_out] = LZc_baseline_multishuffle( ...
                    x_seg, seed + seg*1000 + ch, nBinaryShuffles);

                phase_LZs_values = nan(nPhaseSurrogates, 1);
                phase_raw_values = nan(nPhaseSurrogates, 1);

                for ps = 1:nPhaseSurrogates

                    surr_seed = seed + 1000000 + seg*10000 + ch*100 + ps;
                    x_surr = phase_shuffle_signal(x_seg, surr_seed);

                    [phase_LZs, phase_out] = LZc_baseline_multishuffle( ...
                        x_surr, surr_seed, nBinaryShuffles);

                    phase_LZs_values(ps) = phase_LZs;
                    phase_raw_values(ps) = phase_out.c_orig;

                end

                phase_LZs_mean = mean(phase_LZs_values, 'omitnan');
                phase_LZs_sd = std(phase_LZs_values, 'omitnan');

                if isfinite(phase_LZs_mean) && phase_LZs_mean ~= 0
                    LZsN = real_LZs / phase_LZs_mean;
                else
                    LZsN = NaN;
                end

                row = table( ...
                    string(info.participant), string(info.session), string(info.eyes), string(info.epoch), string(base_name), ...
                    seg, segment_start_sec, segment_end_sec, ...
                    string(chan_labels(ch)), ch, ...
                    real_LZs, real_out.c_orig, real_out.c_shuf_mean, real_out.c_shuf_sd, real_out.TH(1), ...
                    phase_LZs_mean, phase_LZs_sd, mean(phase_raw_values, 'omitnan'), std(phase_raw_values, 'omitnan'), ...
                    LZsN, ...
                    'VariableNames', {'Participant','Session','Eyes','Epoch','File', ...
                    'Segment','SegmentStartSec','SegmentEndSec', ...
                    'Channel','ChannelIndex', ...
                    'LZs','LZs_raw','LZs_binaryShuffleMean','LZs_binaryShuffleSD','Threshold', ...
                    'Phase_LZs_mean','Phase_LZs_sd','Phase_raw_mean','Phase_raw_sd', ...
                    'LZsN'} ...
                );

                LZs_segment_table = [LZs_segment_table; row];

            end
        end

        LZs_segment_table = normalize_lzc_table_ids(LZs_segment_table);

        % Channel summaries across segments.
        channels = unique(LZs_segment_table.Channel, 'stable');

        for ch = 1:numel(channels)

            ch_subset = LZs_segment_table(LZs_segment_table.Channel == channels(ch), :);

            row = table( ...
                string(info.participant), string(info.session), string(info.eyes), string(info.epoch), string(base_name), ...
                channels(ch), ch, ...
                mean(ch_subset.LZs, 'omitnan'), std(ch_subset.LZs, 'omitnan'), ...
                mean(ch_subset.LZsN, 'omitnan'), std(ch_subset.LZsN, 'omitnan'), ...
                mean(ch_subset.LZs_raw, 'omitnan'), std(ch_subset.LZs_raw, 'omitnan'), ...
                mean(ch_subset.Phase_LZs_mean, 'omitnan'), std(ch_subset.Phase_LZs_mean, 'omitnan'), ...
                height(ch_subset), ...
                'VariableNames', {'Participant','Session','Eyes','Epoch','File', ...
                'Channel','ChannelIndex', ...
                'Mean_LZs','SD_LZs','Mean_LZsN','SD_LZsN', ...
                'Mean_LZs_raw','SD_LZs_raw','Mean_Phase_LZs','SD_Phase_LZs', ...
                'N_Segments'} ...
            );

            LZs_channel_summary = [LZs_channel_summary; row];

        end

        LZs_channel_summary = normalize_lzc_table_ids(LZs_channel_summary);

        if save_channel_tables
            writetable(LZs_segment_table, LZs_segment_csv);
        end

        writetable(LZs_channel_summary, LZs_channel_csv);

    end

    %% ---------------- LZc: RANDOM CHANNEL PICKS PER SEGMENT ----------------

    LZc_pick_table = table();
    LZc_segment_summary = table();

    if do_LZc_random_channel_picks

        fprintf('\nComputing LZc using %d random picks of %d channels per segment...\n', ...
            nChannelPicks, nChannelsPerPick);

        if nChannelsPerPick > nChannels
            error('nChannelsPerPick (%d) cannot exceed number of channels (%d).', ...
                nChannelsPerPick, nChannels);
        end

        rng(seed + str2double(regexprep(char(info.participant), '\D', '')));

        % Pre-generate channel picks so every segment uses a reproducible selection set.
        channel_picks = zeros(nChannelPicks, nChannelsPerPick);
        for pick = 1:nChannelPicks
            channel_picks(pick, :) = randperm(nChannels, nChannelsPerPick);
        end

        for seg = 1:nSegments

            idx = seg_starts(seg):(seg_starts(seg) + segment_samples - 1);
            segment_start_sec = (seg_starts(seg) - 1) / srate;
            segment_end_sec = segment_start_sec + segment_sec;

            fprintf('  LZc segment %d/%d\n', seg, nSegments);

            seg_pick_LZc = nan(nChannelPicks, 1);
            seg_pick_LZcN = nan(nChannelPicks, 1);

            for pick = 1:nChannelPicks

                pick_channels = channel_picks(pick, :);
                X_seg_pick = X(pick_channels, idx);

                real_seed = seed + 2000000 + seg*10000 + pick;

                [real_LZc, real_out] = LZc_baseline_multishuffle( ...
                    X_seg_pick, real_seed, nBinaryShuffles);

                phase_LZc_values = nan(nPhaseSurrogates, 1);
                phase_raw_values = nan(nPhaseSurrogates, 1);

                for ps = 1:nPhaseSurrogates

                    surr_seed = seed + 3000000 + seg*100000 + pick*100 + ps;

                    X_surr_pick = zeros(size(X_seg_pick));

                    for pc = 1:size(X_seg_pick, 1)
                        X_surr_pick(pc, :) = phase_shuffle_signal(X_seg_pick(pc, :), surr_seed + pc);
                    end

                    [phase_LZc, phase_out] = LZc_baseline_multishuffle( ...
                        X_surr_pick, surr_seed, nBinaryShuffles);

                    phase_LZc_values(ps) = phase_LZc;
                    phase_raw_values(ps) = phase_out.c_orig;

                end

                phase_LZc_mean = mean(phase_LZc_values, 'omitnan');
                phase_LZc_sd = std(phase_LZc_values, 'omitnan');

                if isfinite(phase_LZc_mean) && phase_LZc_mean ~= 0
                    LZcN = real_LZc / phase_LZc_mean;
                else
                    LZcN = NaN;
                end

                seg_pick_LZc(pick) = real_LZc;
                seg_pick_LZcN(pick) = LZcN;

                channel_pick_labels = strjoin(chan_labels(pick_channels), ",");

                row = table( ...
                    string(info.participant), string(info.session), string(info.eyes), string(info.epoch), string(base_name), ...
                    seg, segment_start_sec, segment_end_sec, ...
                    pick, string(channel_pick_labels), ...
                    real_LZc, real_out.c_orig, real_out.c_shuf_mean, real_out.c_shuf_sd, ...
                    phase_LZc_mean, phase_LZc_sd, mean(phase_raw_values, 'omitnan'), std(phase_raw_values, 'omitnan'), ...
                    LZcN, ...
                    'VariableNames', {'Participant','Session','Eyes','Epoch','File', ...
                    'Segment','SegmentStartSec','SegmentEndSec', ...
                    'Pick','PickedChannels', ...
                    'LZc','LZc_raw','LZc_binaryShuffleMean','LZc_binaryShuffleSD', ...
                    'Phase_LZc_mean','Phase_LZc_sd','Phase_raw_mean','Phase_raw_sd', ...
                    'LZcN'} ...
                );

                LZc_pick_table = [LZc_pick_table; row];

            end

            seg_row = table( ...
                string(info.participant), string(info.session), string(info.eyes), string(info.epoch), string(base_name), ...
                seg, segment_start_sec, segment_end_sec, ...
                mean(seg_pick_LZc, 'omitnan'), std(seg_pick_LZc, 'omitnan'), ...
                mean(seg_pick_LZcN, 'omitnan'), std(seg_pick_LZcN, 'omitnan'), ...
                nChannelPicks, nChannelsPerPick, ...
                'VariableNames', {'Participant','Session','Eyes','Epoch','File', ...
                'Segment','SegmentStartSec','SegmentEndSec', ...
                'Mean_LZc','SD_LZc','Mean_LZcN','SD_LZcN', ...
                'N_ChannelPicks','N_ChannelsPerPick'} ...
            );

            LZc_segment_summary = [LZc_segment_summary; seg_row];

        end

        LZc_pick_table = normalize_lzc_table_ids(LZc_pick_table);
        LZc_segment_summary = normalize_lzc_table_ids(LZc_segment_summary);

        if save_pick_tables
            writetable(LZc_pick_table, LZc_pick_csv);
        end

        if save_segment_tables
            writetable(LZc_segment_summary, LZc_segment_csv);
        end

    end

    %% ---------------- FILE-LEVEL SUMMARY ----------------

    file_summary = table( ...
        string(info.participant), string(info.session), string(info.eyes), string(info.epoch), string(base_name), ...
        nChannels, nSamples, srate, segment_sec, nSegments, nBinaryShuffles, nPhaseSurrogates, ...
        'VariableNames', {'Participant','Session','Eyes','Epoch','File', ...
        'N_Channels','N_Samples','SamplingRate','SegmentSec','N_Segments', ...
        'N_BinaryShuffles','N_PhaseSurrogates'} ...
    );

    if do_LZs_by_channel && ~isempty(LZs_segment_table)
        file_summary.Mean_LZs = mean(LZs_segment_table.LZs, 'omitnan');
        file_summary.SD_LZs = std(LZs_segment_table.LZs, 'omitnan');
        file_summary.Mean_LZsN = mean(LZs_segment_table.LZsN, 'omitnan');
        file_summary.SD_LZsN = std(LZs_segment_table.LZsN, 'omitnan');
        file_summary.Mean_LZs_raw = mean(LZs_segment_table.LZs_raw, 'omitnan');
        file_summary.Mean_Phase_LZs = mean(LZs_segment_table.Phase_LZs_mean, 'omitnan');
    else
        file_summary.Mean_LZs = NaN;
        file_summary.SD_LZs = NaN;
        file_summary.Mean_LZsN = NaN;
        file_summary.SD_LZsN = NaN;
        file_summary.Mean_LZs_raw = NaN;
        file_summary.Mean_Phase_LZs = NaN;
    end

    if do_LZc_random_channel_picks && ~isempty(LZc_segment_summary)
        file_summary.Mean_LZc = mean(LZc_segment_summary.Mean_LZc, 'omitnan');
        file_summary.SD_LZc = std(LZc_segment_summary.Mean_LZc, 'omitnan');
        file_summary.Mean_LZcN = mean(LZc_segment_summary.Mean_LZcN, 'omitnan');
        file_summary.SD_LZcN = std(LZc_segment_summary.Mean_LZcN, 'omitnan');
        file_summary.N_ChannelPicks = nChannelPicks;
        file_summary.N_ChannelsPerPick = nChannelsPerPick;
    else
        file_summary.Mean_LZc = NaN;
        file_summary.SD_LZc = NaN;
        file_summary.Mean_LZcN = NaN;
        file_summary.SD_LZcN = NaN;
        file_summary.N_ChannelPicks = NaN;
        file_summary.N_ChannelsPerPick = NaN;
    end

    file_summary = normalize_lzc_table_ids(file_summary);

    writetable(file_summary, summary_csv);

    save(mat_path, ...
        'file_summary', ...
        'LZs_channel_summary', ...
        'LZc_segment_summary', ...
        'info', ...
        'set_path', ...
        'segment_sec', ...
        'segment_overlap_sec', ...
        'nBinaryShuffles', ...
        'nPhaseSurrogates', ...
        'nChannelPicks', ...
        'nChannelsPerPick', ...
        'seed', ...
        '-v7.3');

    all_file_summary = [all_file_summary; file_summary];
    all_LZs_channel_summary = [all_LZs_channel_summary; LZs_channel_summary];
    all_LZc_segment_summary = [all_LZc_segment_summary; LZc_segment_summary];

    fprintf('\nSaved file summary:\n%s\n', summary_csv);

end

%% ---------------- COMBINED OUTPUTS ----------------

all_file_summary = normalize_lzc_table_ids(all_file_summary);
all_LZs_channel_summary = normalize_lzc_table_ids(all_LZs_channel_summary);
all_LZc_segment_summary = normalize_lzc_table_ids(all_LZc_segment_summary);

combined_file_csv = fullfile(output_root, 'ALL_paper_method_file_summary.csv');
combined_LZs_channel_csv = fullfile(output_root, 'ALL_paper_method_LZs_channel_summary.csv');
combined_LZc_segment_csv = fullfile(output_root, 'ALL_paper_method_LZc_segment_summary.csv');
combined_mat = fullfile(output_root, 'ALL_paper_method_combined_outputs.mat');

writetable(all_file_summary, combined_file_csv);

if ~isempty(all_LZs_channel_summary)
    writetable(all_LZs_channel_summary, combined_LZs_channel_csv);
end

if ~isempty(all_LZc_segment_summary)
    writetable(all_LZc_segment_summary, combined_LZc_segment_csv);
end

save(combined_mat, ...
    'all_file_summary', ...
    'all_LZs_channel_summary', ...
    'all_LZc_segment_summary', ...
    'segment_sec', ...
    'segment_overlap_sec', ...
    'nBinaryShuffles', ...
    'nPhaseSurrogates', ...
    'nChannelPicks', ...
    'nChannelsPerPick', ...
    'seed', ...
    '-v7.3');

fprintf('\nSaved combined file summary:\n%s\n', combined_file_csv);
fprintf('Saved combined MAT:\n%s\n', combined_mat);
fprintf('\nPaper-method LZc/LZcN pipeline complete.\n');
