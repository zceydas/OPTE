%% lzc_heatmap.m
clear; clc;
format long g

%% ---------------- EEGLAB SETUP ----------------

eeglab_path = '/Users/luisluna/Documents/MATLAB/eeglab2025.0.0';

if exist(eeglab_path, 'dir')
    addpath(eeglab_path);
    eeglab;
else
    error('EEGLAB path not found: %s', eeglab_path);
end

%% ---------------- USER SETTINGS ----------------

seed = 0;

nShuffles = 10;
plot_lzc_column = "LZc_10";

lzc_color_limits = [0.45 0.85];
diff_color_limits = [-0.15 0.15];

do_spatial_lzc = true;
do_temporal_lzc = true;
do_eo_minus_ec = true;
do_temporal_eo_minus_ec = true;
do_participant_dosing_summary = true;
do_group_averages = true;

temporal_window_sec = 10;
temporal_step_sec = 2;
temporal_nShuffles = 10;
temporal_plot_lzc_column = "LZc_10";

session_order = {'baseline', 'dosing', '1week', '2week', '1month'};
session_labels = {'Baseline', 'Dosing', '1 Week', '2 Week', '1 Month'};

eyes_order = {'EO', 'EC'};
eyes_labels = {'Eyes Open', 'Eyes Closed'};

%% ---------------- SELECT INPUT FOLDER ----------------

input_root = uigetdir(pwd, 'Select folder containing EEGLAB .set files');

if isequal(input_root, 0)
    disp('Folder selection cancelled.');
    return;
end

output_root = fullfile(input_root, 'LZc_Results');

if ~exist(output_root, 'dir')
    mkdir(output_root);
end

%% ---------------- FIND FILES ----------------

set_files = dir(fullfile(input_root, '**', 'PostICA_*.set'));

if isempty(set_files)
    error('No PostICA_*.set files found under: %s', input_root);
end

fprintf('\nFound %d PostICA .set files.\n', numel(set_files));

all_results = table();

template_chanlocs = [];
template_labels = strings(0,1);

%% ---------------- PROCESS EACH FILE ----------------

for f = 1:numel(set_files)

    file = set_files(f).name;
    path = set_files(f).folder;
    set_path = fullfile(path, file);
    [~, base_name, ~] = fileparts(file);

    fprintf('\n====================================================\n');
    fprintf('Processing %d/%d:\n%s\n', f, numel(set_files), set_path);

    info = parse_lzc_filename(base_name);

    if info.skip
        fprintf('Skipping file because participant/session/eyes condition could not be parsed.\n');
        continue;
    end

    info.participant = normalize_participant_id(info.participant);

    participant_dir = fullfile(output_root, char(info.participant));

    if ~exist(participant_dir, 'dir')
        mkdir(participant_dir);
    end

    file_tag = sprintf('%s_%s_%s_%s_%s', ...
        char(info.participant), char(info.session), char(info.eyes), char(info.epoch), base_name);

    csv_path = fullfile(participant_dir, [file_tag '_LZc_by_channel.csv']);
    mat_path = fullfile(participant_dir, [file_tag '_LZc_outputs.mat']);

    channel_heatmap_path = fullfile(participant_dir, ...
        [file_tag '_' char(plot_lzc_column) '_channel_heatmap.png']);

    topoplot_path = fullfile(participant_dir, ...
        [file_tag '_' char(plot_lzc_column) '_topoplot.png']);

    temporal_dir = fullfile(participant_dir, 'Temporal_LZc');

    if ~exist(temporal_dir, 'dir')
        mkdir(temporal_dir);
    end

    temporal_mat_path = fullfile(temporal_dir, ...
        [file_tag '_' char(temporal_plot_lzc_column) '_temporal_LZc.mat']);

    temporal_fig_path = fullfile(temporal_dir, ...
        [file_tag '_' char(temporal_plot_lzc_column) '_temporal_LZc_heatmap.png']);

    need_load_EEG = false;

    if do_spatial_lzc
        if ~exist(csv_path, 'file') || ~exist(mat_path, 'file') || ...
                ~exist(channel_heatmap_path, 'file') || ~exist(topoplot_path, 'file')
            need_load_EEG = true;
        end
    end

    if do_temporal_lzc
        if ~exist(temporal_mat_path, 'file') || ~exist(temporal_fig_path, 'file')
            need_load_EEG = true;
        end
    end

    if need_load_EEG
        EEG = pop_loadset('filename', file, 'filepath', path);
        X = double(EEG.data);
        [nChannels, ~] = size(X);
        chan_labels = get_channel_labels(EEG, nChannels);

        if isempty(template_chanlocs) && isfield(EEG, 'chanlocs') && numel(EEG.chanlocs) >= nChannels
            template_chanlocs = EEG.chanlocs;
            template_labels = chan_labels;
        end
    else
        EEG = [];
        X = [];
        nChannels = [];
        chan_labels = [];
        fprintf('All file-level outputs already exist. Loading CSV only.\n');
    end

    %% ---------------- SPATIAL LZc ----------------

    if do_spatial_lzc

        if exist(csv_path, 'file') && exist(mat_path, 'file')
            fprintf('Spatial LZc already computed. Loading CSV:\n%s\n', csv_path);
            results_table = readtable(csv_path, 'TextType', 'string');
            results_table = normalize_lzc_table_ids(results_table);

        else
            fprintf('Computing spatial channel-wise LZc with %d shuffles...\n', nShuffles);

            LZc_10 = nan(nChannels, 1);
            c_orig_values = nan(nChannels, 1);
            c_shuf_mean_10 = nan(nChannels, 1);
            c_shuf_sd_10 = nan(nChannels, 1);
            threshold_values = nan(nChannels, 1);
            all_out = cell(nChannels, 1);

            for ch = 1:nChannels

                fprintf('  Channel %d/%d: %s\n', ch, nChannels, chan_labels(ch));

                X_ch = X(ch, :);

                [lz_val, out] = LZc_baseline_multishuffle(X_ch, seed + ch, nShuffles);

                LZc_10(ch) = lz_val;
                c_orig_values(ch) = out.c_orig;
                c_shuf_mean_10(ch) = out.c_shuf_mean;
                c_shuf_sd_10(ch) = out.c_shuf_sd;
                threshold_values(ch) = out.TH(1);
                all_out{ch} = out;

            end

            results_table = table( ...
                repmat(string(info.participant), nChannels, 1), ...
                repmat(string(info.session), nChannels, 1), ...
                repmat(string(info.eyes), nChannels, 1), ...
                repmat(string(info.epoch), nChannels, 1), ...
                repmat(string(base_name), nChannels, 1), ...
                chan_labels, ...
                LZc_10, ...
                c_orig_values, ...
                c_shuf_mean_10, ...
                c_shuf_sd_10, ...
                threshold_values, ...
                'VariableNames', {'Participant','Session','Eyes','Epoch','File','Channel', ...
                'LZc_10','c_orig','c_shuf_mean_10','c_shuf_sd_10','Threshold'} ...
            );

            results_table = normalize_lzc_table_ids(results_table);

            writetable(results_table, csv_path);

            save(mat_path, ...
                'results_table', ...
                'LZc_10', ...
                'c_orig_values', ...
                'c_shuf_mean_10', ...
                'c_shuf_sd_10', ...
                'threshold_values', ...
                'all_out', ...
                'chan_labels', ...
                'set_path', ...
                'seed', ...
                'nShuffles', ...
                'lzc_color_limits', ...
                '-v7.3');

            fprintf('Saved spatial LZc CSV:\n%s\n', csv_path);
            fprintf('Saved spatial LZc MAT:\n%s\n', mat_path);

        end

        if ~exist(channel_heatmap_path, 'file')
            fprintf('Creating channel heatmap...\n');

            plot_LZc_values = results_table.(plot_lzc_column);

            plot_channel_heatmap(plot_LZc_values, results_table.Channel, ...
                [file_tag ' Channel-wise LZc using ' char(plot_lzc_column)], ...
                channel_heatmap_path, ...
                lzc_color_limits);
        else
            fprintf('Channel heatmap already exists. Skipping:\n%s\n', channel_heatmap_path);
        end

        if ~exist(topoplot_path, 'file')
            fprintf('Creating spatial topoplot...\n');

            if isempty(EEG)
                EEG = pop_loadset('filename', file, 'filepath', path);
            end

            plot_LZc_values = results_table.(plot_lzc_column);

            if exist('topoplot', 'file') == 2 && isfield(EEG, 'chanlocs')
                plot_lzc_topoplot(plot_LZc_values, EEG.chanlocs, ...
                    [file_tag ' LZc Topoplot using ' char(plot_lzc_column)], ...
                    topoplot_path, ...
                    lzc_color_limits);
            else
                fprintf('Skipping topoplot because topoplot or EEG.chanlocs is unavailable.\n');
            end
        else
            fprintf('Topoplot already exists. Skipping:\n%s\n', topoplot_path);
        end

        all_results = [all_results; results_table];

    end

    %% ---------------- TEMPORAL LZc ----------------

    if do_temporal_lzc

        if exist(temporal_mat_path, 'file') && exist(temporal_fig_path, 'file')
            fprintf('Temporal LZc outputs already exist. Skipping temporal step:\n%s\n', temporal_mat_path);

        else
            fprintf('Computing temporal LZc for %s...\n', file_tag);

            if isempty(EEG)
                EEG = pop_loadset('filename', file, 'filepath', path);
                X = double(EEG.data);
                [nChannels, ~] = size(X);
                chan_labels = get_channel_labels(EEG, nChannels);
            end

            [temporal_LZc, temporal_times_sec] = temporal_lzc_by_channel( ...
                X, EEG.srate, temporal_window_sec, temporal_step_sec, seed, temporal_nShuffles);

            save(temporal_mat_path, ...
                'temporal_LZc', ...
                'temporal_times_sec', ...
                'chan_labels', ...
                'temporal_window_sec', ...
                'temporal_step_sec', ...
                'temporal_nShuffles', ...
                'temporal_plot_lzc_column', ...
                'set_path', ...
                '-v7.3');

            temporal_lzc_heatmap( ...
                temporal_LZc, temporal_times_sec, chan_labels, ...
                [file_tag ' Temporal LZc using ' char(temporal_plot_lzc_column)], ...
                temporal_fig_path, ...
                lzc_color_limits);

            fprintf('Saved temporal LZc MAT:\n%s\n', temporal_mat_path);
            fprintf('Saved temporal LZc heatmap:\n%s\n', temporal_fig_path);

        end

    end

end

%% ---------------- SAVE COMBINED TABLE ----------------

all_results = normalize_lzc_table_ids(all_results);

if isempty(all_results)
    error('No valid spatial LZc results were available. Check parsing and file outputs.');
end

combined_csv = fullfile(output_root, 'ALL_LZc_by_channel.csv');
combined_mat = fullfile(output_root, 'ALL_LZc_by_channel.mat');

writetable(all_results, combined_csv);

save(combined_mat, ...
    'all_results', ...
    'seed', ...
    'nShuffles', ...
    'plot_lzc_column', ...
    'lzc_color_limits', ...
    '-v7.3');

fprintf('\nSaved combined table:\n%s\n', combined_csv);

%% ---------------- ENSURE TEMPLATE CHANNEL LOCATIONS EXIST ----------------

if isempty(template_chanlocs) || isempty(template_labels)

    fprintf('\nTemplate channel locations not loaded yet. Loading first valid .set file for template...\n');

    template_loaded = false;

    for tf = 1:numel(set_files)

        temp_file = set_files(tf).name;
        temp_path = set_files(tf).folder;

        temp_EEG = pop_loadset('filename', temp_file, 'filepath', temp_path);

        if isfield(temp_EEG, 'chanlocs') && ~isempty(temp_EEG.chanlocs)

            template_chanlocs = temp_EEG.chanlocs;
            template_labels = get_channel_labels(temp_EEG, temp_EEG.nbchan);

            template_loaded = true;

            fprintf('Loaded template channel locations from:\n%s\n', fullfile(temp_path, temp_file));
            break;
        end
    end

    if ~template_loaded
        warning('Could not load template channel locations. Topoplot-based summaries may fail or be skipped.');
    end
end

%% ---------------- SPATIAL EO MINUS EC DIFFERENCE MAPS ----------------

if do_eo_minus_ec

    diff_dir = fullfile(output_root, 'EO_minus_EC_Differences');

    if ~exist(diff_dir, 'dir')
        mkdir(diff_dir);
    end

    diff_csv = fullfile(diff_dir, 'ALL_EO_minus_EC_LZc_by_channel.csv');
    diff_mat = fullfile(diff_dir, 'ALL_EO_minus_EC_LZc_by_channel.mat');

    if exist(diff_csv, 'file') && exist(diff_mat, 'file')
        fprintf('\nSpatial EO minus EC difference table already exists. Loading:\n%s\n', diff_csv);
        diff_results = readtable(diff_csv, 'TextType', 'string');
        diff_results = normalize_lzc_table_ids(diff_results);
    else
        fprintf('\nComputing spatial EO minus EC differences...\n');

        diff_results = eo_minus_ec_differences(all_results, plot_lzc_column);
        diff_results = normalize_lzc_table_ids(diff_results);

        writetable(diff_results, diff_csv);

        save(diff_mat, ...
            'diff_results', ...
            'plot_lzc_column', ...
            'diff_color_limits', ...
            '-v7.3');

        fprintf('Saved spatial EO minus EC difference table:\n%s\n', diff_csv);
    end

    participants = unique(string(diff_results.Participant), 'stable');

    for p = 1:numel(participants)

        participant_id = normalize_participant_id(participants(p));
        participant_id_clean = char(participant_id);

        participant_diff_dir = fullfile(diff_dir, participant_id_clean);

        if ~exist(participant_diff_dir, 'dir')
            mkdir(participant_diff_dir);
        end

        participant_diff_path = fullfile(participant_diff_dir, ...
            ['Participant_' participant_id_clean '_EO_minus_EC_' char(plot_lzc_column) '.png']);

        if exist(participant_diff_path, 'file')
            fprintf('Participant spatial EO minus EC map already exists. Skipping:\n%s\n', participant_diff_path);
        else
            plot_participant_difference_maps( ...
                diff_results, participant_id, ...
                template_labels, template_chanlocs, ...
                session_order, session_labels, ...
                plot_lzc_column, diff_color_limits, ...
                participant_diff_dir);
        end

    end

end

%% ---------------- EO MINUS EC BAR + CHANNEL SUMMARY ----------------

summary_bar_dir = fullfile(output_root, 'EO_minus_EC_Group_Summary');

summary_bar_path = fullfile(summary_bar_dir, ...
    ['EO_minus_EC_group_bar_summary_' char(plot_lzc_column) '.png']);

summary_heatmap_path = fullfile(summary_bar_dir, ...
    ['EO_minus_EC_channel_by_condition_heatmap_' char(plot_lzc_column) '.png']);

if exist(summary_bar_path, 'file') && exist(summary_heatmap_path, 'file')
    fprintf('EO minus EC group summary figures already exist. Skipping.\n');
else
    plot_eo_minus_ec_group_bar_and_channel_summary( ...
        diff_results, ...
        plot_lzc_column, ...
        summary_bar_dir);
end

%% ---------------- TEMPORAL EO MINUS EC DIFFERENCES ----------------

if do_temporal_eo_minus_ec

    temporal_diff_dir = fullfile(output_root, 'Temporal_EO_minus_EC_Differences');

    if ~exist(temporal_diff_dir, 'dir')
        mkdir(temporal_diff_dir);
    end

    fprintf('\nComputing participant-level temporal EO minus EC differences...\n');

    temporal_pair_index = temporal_eo_minus_ec_pair_differences( ...
        all_results, ...
        output_root, ...
        temporal_plot_lzc_column, ...
        diff_color_limits, ...
        temporal_diff_dir);

    pair_index_csv = fullfile(temporal_diff_dir, 'TEMPORAL_EO_minus_EC_pair_index.csv');
    writetable(temporal_pair_index, pair_index_csv);

    pair_index_mat = fullfile(temporal_diff_dir, 'TEMPORAL_EO_minus_EC_pair_index.mat');
    save(pair_index_mat, 'temporal_pair_index', '-v7.3');

    fprintf('\nComputing group/session average temporal EO minus EC differences...\n');

    temporal_group_dir = fullfile(temporal_diff_dir, 'Group_Averages');

    if ~exist(temporal_group_dir, 'dir')
        mkdir(temporal_group_dir);
    end

    temporal_group_index = temporal_eo_minus_ec_group_averages( ...
        temporal_pair_index, ...
        session_order, ...
        session_labels, ...
        temporal_plot_lzc_column, ...
        diff_color_limits, ...
        temporal_group_dir);

    group_index_csv = fullfile(temporal_group_dir, 'TEMPORAL_EO_minus_EC_group_index.csv');
    writetable(temporal_group_index, group_index_csv);

    group_index_mat = fullfile(temporal_group_dir, 'TEMPORAL_EO_minus_EC_group_index.mat');
    save(group_index_mat, 'temporal_group_index', '-v7.3');

end

%% ---------------- PARTICIPANT DOSING 2 x 4 SUMMARY FIGURES ----------------

if do_participant_dosing_summary

    participant_summary_dir = fullfile(output_root, 'Participant_Summaries');

    if ~exist(participant_summary_dir, 'dir')
        mkdir(participant_summary_dir);
    end

    participants = unique(string(all_results.Participant), 'stable');

    for p = 1:numel(participants)

        participant_id = normalize_participant_id(participants(p));
        participant_id_clean = char(participant_id);

        dosing_summary_path = fullfile(participant_summary_dir, ...
            ['Participant_' participant_id_clean '_dosing_2x4_' char(plot_lzc_column) '.png']);

        if exist(dosing_summary_path, 'file')
            fprintf('Participant dosing 2x4 summary already exists. Skipping:\n%s\n', dosing_summary_path);
        else
            plot_participant_dosing_2x4( ...
                all_results, participant_id, ...
                template_labels, template_chanlocs, ...
                eyes_order, eyes_labels, ...
                plot_lzc_column, lzc_color_limits, ...
                participant_summary_dir);
        end

    end

end

%% ---------------- GROUP AVERAGES ----------------

if do_group_averages

    avg_dir = fullfile(output_root, 'Group_Averages');

    if ~exist(avg_dir, 'dir')
        mkdir(avg_dir);
    end

    avg_csv = fullfile(avg_dir, 'ALL_GROUP_AVERAGE_LZc.csv');
    avg_mat = fullfile(avg_dir, 'ALL_GROUP_AVERAGE_LZc.mat');

    fprintf('\nComputing group averages...\n');

    avg_results = table();

    for s = 1:numel(session_order)

        session_name = session_order{s};

        for e = 1:numel(eyes_order)

            eyes_name = eyes_order{e};

            subset = all_results( ...
                all_results.Session == string(session_name) & ...
                all_results.Eyes == string(eyes_name), :);

            if isempty(subset)
                fprintf('\nNo data for %s %s. Skipping group average.\n', session_name, eyes_name);
                continue;
            end

            channels = unique(subset.Channel, 'stable');
            nChannels = numel(channels);

            Mean_LZc_10 = nan(nChannels, 1);
            SD_LZc_10 = nan(nChannels, 1);
            n_participants = nan(nChannels, 1);

            for ch = 1:nChannels

                ch_name = channels(ch);
                ch_subset = subset(subset.Channel == ch_name, :);

                participant_means = groupsummary(ch_subset, 'Participant', 'mean', 'LZc_10');

                Mean_LZc_10(ch) = mean(participant_means.mean_LZc_10, 'omitnan');
                SD_LZc_10(ch) = std(participant_means.mean_LZc_10, 'omitnan');
                n_participants(ch) = height(participant_means);

            end

            avg_table = table( ...
                repmat(string(session_name), nChannels, 1), ...
                repmat(string(eyes_name), nChannels, 1), ...
                channels, ...
                Mean_LZc_10, ...
                SD_LZc_10, ...
                n_participants, ...
                'VariableNames', {'Session','Eyes','Channel', ...
                'Mean_LZc_10','SD_LZc_10','N_Participants'} ...
            );

            avg_table = normalize_lzc_table_ids(avg_table);
            avg_results = [avg_results; avg_table];

            tag = sprintf('GROUP_%s_%s', session_name, eyes_name);

            group_csv_path = fullfile(avg_dir, [tag '_average_LZc.csv']);
            group_heatmap_path = fullfile(avg_dir, ...
                [tag '_' char(plot_lzc_column) '_average_channel_heatmap.png']);
            group_topoplot_path = fullfile(avg_dir, ...
                [tag '_' char(plot_lzc_column) '_average_topoplot.png']);

            writetable(avg_table, group_csv_path);

            if ~exist(group_heatmap_path, 'file')
                avg_plot_values = avg_table.Mean_LZc_10;

                plot_channel_heatmap(avg_plot_values, channels, ...
                    [tag ' Average Channel-wise LZc using ' char(plot_lzc_column)], ...
                    group_heatmap_path, ...
                    lzc_color_limits);
            end

            if ~exist(group_topoplot_path, 'file') && exist('topoplot', 'file') == 2 && ~isempty(template_chanlocs)

                [aligned_values, aligned_chanlocs] = align_average_to_template( ...
                    avg_table, template_labels, template_chanlocs, 'Mean_LZc_10');

                plot_lzc_topoplot(aligned_values, aligned_chanlocs, ...
                    [tag ' Average LZc Topoplot using ' char(plot_lzc_column)], ...
                    group_topoplot_path, ...
                    lzc_color_limits);
            end

        end

    end

    writetable(avg_results, avg_csv);

    save(avg_mat, ...
        'avg_results', ...
        'plot_lzc_column', ...
        'lzc_color_limits', ...
        '-v7.3');

    fprintf('\nSaved group averages:\n%s\n', avg_csv);

    combined_heatmap_path = fullfile(avg_dir, ...
        ['GROUP_ALL_SESSIONS_2x5_' char(plot_lzc_column) '_channel_heatmaps.png']);

    if ~exist(combined_heatmap_path, 'file')
        plot_group_2x5_channel_heatmaps( ...
            avg_results, session_order, session_labels, eyes_order, eyes_labels, ...
            plot_lzc_column, combined_heatmap_path, lzc_color_limits);
    else
        fprintf('Combined 2x5 channel heatmap already exists. Skipping:\n%s\n', combined_heatmap_path);
    end

    combined_topoplot_path = fullfile(avg_dir, ...
        ['GROUP_ALL_SESSIONS_2x5_' char(plot_lzc_column) '_topoplots.png']);

    if ~exist(combined_topoplot_path, 'file') && exist('topoplot', 'file') == 2 && ~isempty(template_chanlocs)
        plot_group_2x5_topoplots( ...
            avg_results, template_labels, template_chanlocs, ...
            session_order, session_labels, eyes_order, eyes_labels, ...
            plot_lzc_column, combined_topoplot_path, lzc_color_limits);
    elseif exist(combined_topoplot_path, 'file')
        fprintf('Combined 2x5 topoplot already exists. Skipping:\n%s\n', combined_topoplot_path);
    else
        fprintf('Skipping combined 2x5 topoplot: topoplot or template channel locations unavailable.\n');
    end

end

fprintf('\nAll LZc batch processing complete.\n');