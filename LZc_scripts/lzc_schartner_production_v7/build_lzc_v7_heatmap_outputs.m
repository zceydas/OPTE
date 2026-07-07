%% build_lzc_v7_heatmap_outputs.m
% Build all V7-derived outputs:
%   1) channel-wise heatmap/topoplot outputs for LZs and LZsN
%   2) all-channel file-level outputs for LZc and LZcN
%
% Required inputs:
%   LZc_Schartner_ProductionV7_Results/CSV/PUBLIC_long_LZs_by_channel_window.csv
%   LZc_Schartner_ProductionV7_Results/CSV/PUBLIC_long_LZc_by_window.csv
%   LZc_Schartner_ProductionV7_Results/CSV/PUBLIC_full_recording_LZs_by_channel.csv
%   LZc_Schartner_ProductionV7_Results/CSV/PUBLIC_full_recording_LZc_by_file.csv
%
% Run run_lzc_v7_full_recording_supplement.m first if the two full-recording
% CSVs do not exist yet.

clear; clc;
format long g

%% ---------------- PATH CONFIGURATION ----------------
% Set USE_HARDCODED_PATHS=true for reproducible batch runs.
% Leave it false to choose the V7 results folder interactively.
USE_HARDCODED_PATHS = false;

% Optional hardcoded paths. Edit these for your machine if needed.
HARD_CODED_EEGLAB_PATH = '';          % e.g., '/path/to/eeglab2026.0.0'
HARD_CODED_V7_RESULTS_ROOT = '';      % V7 results folder, its parent input folder, or its CSV folder
HARD_CODED_LEGACY_LZC_HELPER_DIR = ''; % optional helper folder if your setup keeps helpers elsewhere

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir, '-begin');

eeglab_path = char(HARD_CODED_EEGLAB_PATH);
legacy_lzc_scripts_dir = char(HARD_CODED_LEGACY_LZC_HELPER_DIR);

%% ---------------- SETTINGS ----------------

data_sources = ["window10s", "window2s", "full"];
channel_metrics = ["LZs", "LZsN"];
all_channel_metrics = ["LZc", "LZcN"];

session_order = {'baseline', 'dosing', '1week', '2week', '1month'};
session_labels = {'Baseline', 'Dosing', '1 Week', '2 Week', '1 Month'};

eyes_order = {'EO', 'EC'};
eyes_labels = {'Eyes Open', 'Eyes Closed'};

%% ---------------- SETUP ----------------

add_optional_path(eeglab_path, 'EEGLAB');
add_optional_path(legacy_lzc_scripts_dir, 'legacy LZc helper folder');

eeglab;

old_fig_visibility = get(groot, 'defaultFigureVisible');
set(groot, 'defaultFigureVisible', 'off');
cleanup_obj = onCleanup(@() set(groot, 'defaultFigureVisible', old_fig_visibility));

if USE_HARDCODED_PATHS
    selected_root = char(HARD_CODED_V7_RESULTS_ROOT);
    if isempty(selected_root)
        error('HARD_CODED_V7_RESULTS_ROOT is empty. Set it or use interactive folder selection.');
    end
else
    selected_root = uigetdir(pwd, ['Select folder containing V7 results, ' ...
        'LZc_Schartner_ProductionV7_Results, or its CSV folder']);
    if isequal(selected_root, 0)
        disp('Folder selection cancelled.');
        return;
    end
end

if ~exist(selected_root, 'dir')
    error('Selected folder not found: %s', selected_root);
end

[input_root, v7_root, csv_root] = resolve_v7_paths(selected_root);

if ~exist(csv_root, 'dir')
    error('V7 CSV folder not found: %s', csv_root);
end

fprintf('\nUsing input root:\n%s\n', input_root);
fprintf('Using V7 root:\n%s\n', v7_root);
fprintf('Using CSV root:\n%s\n', csv_root);

[template_chanlocs, template_labels] = load_template_chanlocs(input_root);

%% ---------------- CHANNEL-WISE LZs/LZsN OUTPUTS ----------------

for ds = 1:numel(data_sources)
    for m = 1:numel(channel_metrics)

        data_source = data_sources(ds);
        metric = channel_metrics(m);

        fprintf('\n====================================================\n');
        fprintf('Building channel-wise outputs for %s / %s\n', char(data_source), char(metric));

        channel_table = make_channel_table(csv_root, data_source);
        channel_table = normalize_lzc_table_ids(channel_table);

        if isempty(channel_table)
            error('No channel rows available for %s.', char(data_source));
        end

        metric_values = channel_table.(char(metric));
        metric_color_limits = percentile_limits(metric_values, [2.5 97.5]);

        out_root = fullfile(v7_root, 'V7_Derived_Heatmaps', char(data_source), char(metric));
        out_csv = fullfile(out_root, 'CSV');
        out_fig = fullfile(out_root, 'Figures');
        per_file_fig_dir = fullfile(out_fig, 'Per_File');
        group_fig_dir = fullfile(out_fig, 'Group_Averages');
        diff_dir = fullfile(out_root, 'EO_minus_EC_Differences');

        mkdir_if_needed(out_csv);
        mkdir_if_needed(per_file_fig_dir);
        mkdir_if_needed(group_fig_dir);
        mkdir_if_needed(diff_dir);

        writetable(channel_table, fullfile(out_csv, ...
            ['ALL_' char(data_source) '_' char(metric) '_by_channel.csv']));

        build_per_file_channel_figures(channel_table, data_source, metric, metric_color_limits, ...
            per_file_fig_dir, template_labels, template_chanlocs);

        diff_table = make_eo_minus_ec_channel_table(channel_table, metric);
        writetable(diff_table, fullfile(diff_dir, ...
            ['ALL_EO_minus_EC_' char(data_source) '_' char(metric) '_by_channel.csv']));

        group_table = make_group_average_table(channel_table, metric, session_order, eyes_order);
        writetable(group_table, fullfile(out_csv, ...
            ['ALL_GROUP_AVERAGE_' char(data_source) '_' char(metric) '.csv']));

        build_group_channel_figures(group_table, data_source, metric, metric_color_limits, ...
            group_fig_dir, template_labels, template_chanlocs, ...
            session_order, session_labels, eyes_order, eyes_labels);

        fprintf('Finished channel-wise %s / %s:\n%s\n', char(data_source), char(metric), out_root);

    end
end

%% ---------------- ALL-CHANNEL LZc/LZcN OUTPUTS ----------------

for ds = 1:numel(data_sources)
    for m = 1:numel(all_channel_metrics)

        data_source = data_sources(ds);
        metric = all_channel_metrics(m);

        fprintf('\n====================================================\n');
        fprintf('Building all-channel outputs for %s / %s\n', char(data_source), char(metric));

        file_table = make_lzc_file_table(csv_root, data_source);
        file_table = normalize_lzc_table_ids(file_table);

        if isempty(file_table)
            error('No all-channel rows available for %s.', char(data_source));
        end

        out_root = fullfile(v7_root, 'V7_Derived_AllChannel_LZc', char(data_source), char(metric));
        out_csv = fullfile(out_root, 'CSV');
        out_fig = fullfile(out_root, 'Figures');
        diff_dir = fullfile(out_root, 'EO_minus_EC_Differences');

        mkdir_if_needed(out_csv);
        mkdir_if_needed(out_fig);
        mkdir_if_needed(diff_dir);

        writetable(file_table, fullfile(out_csv, ...
            ['ALL_' char(data_source) '_' char(metric) '_by_file.csv']));

        diff_table = make_eo_minus_ec_file_table(file_table, metric);
        writetable(diff_table, fullfile(diff_dir, ...
            ['ALL_EO_minus_EC_' char(data_source) '_' char(metric) '_by_file.csv']));

        group_table = make_file_group_summary(diff_table, metric, session_order);
        writetable(group_table, fullfile(out_csv, ...
            ['GROUP_EO_minus_EC_' char(data_source) '_' char(metric) '_summary.csv']));

        plot_file_group_summary(group_table, data_source, metric, out_fig);

        fprintf('Finished all-channel %s / %s:\n%s\n', char(data_source), char(metric), out_root);

    end
end

fprintf('\nAll V7-derived LZs/LZsN and LZc/LZcN outputs complete.\n');

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

function [input_root, v7_root, csv_root] = resolve_v7_paths(selected_root)
selected_root = char(selected_root);
[parent_dir, folder_name] = fileparts(selected_root);

if strcmp(folder_name, 'CSV')
    csv_root = selected_root;
    v7_root = parent_dir;
    input_root = fileparts(v7_root);
elseif strcmp(folder_name, 'LZc_Schartner_ProductionV7_Results')
    v7_root = selected_root;
    csv_root = fullfile(v7_root, 'CSV');
    input_root = fileparts(v7_root);
else
    input_root = selected_root;
    v7_root = fullfile(input_root, 'LZc_Schartner_ProductionV7_Results');
    csv_root = fullfile(v7_root, 'CSV');
end
end

function T = make_channel_table(csv_root, data_source)
switch string(data_source)
    case "full"
        p = fullfile(csv_root, 'PUBLIC_full_recording_LZs_by_channel.csv');
        if ~exist(p, 'file')
            error(['Missing full-recording channel CSV:\n%s\n\n' ...
                'Run run_lzc_v7_full_recording_supplement.m first.'], p);
        end
        T = readtable(p, 'TextType', 'string');
        T.DataSource = repmat("full", height(T), 1);

    case {"window10s", "window2s"}
        p = fullfile(csv_root, 'PUBLIC_long_LZs_by_channel_window.csv');
        if ~exist(p, 'file')
            error('V7 public channel-window CSV not found: %s', p);
        end
        W = readtable(p, 'TextType', 'string');

        if string(data_source) == "window10s"
            window_sec = 10;
        else
            window_sec = 2;
        end

        W = W(asnum_local(W.WindowLengthSec) == window_sec, :);
        group_vars = {'Participant','Session','Eyes','Epoch','File','Channel','ChannelIndex'};
        T = groupsummary(W, group_vars, 'mean', ...
            {'LZs','LZsN','RawLZs','BinaryShuffleMeanRawLZs','PhaseRawLZsMean', ...
            'Threshold','PropOnes','NTransitions','StringLength'});
        T.GroupCount = [];
        T.Properties.VariableNames = erase(T.Properties.VariableNames, 'mean_');
        T.WindowLengthSec = repmat(double(window_sec), height(T), 1);
        T.DataSource = repmat(string(data_source), height(T), 1);
        T = movevars(T, {'WindowLengthSec','DataSource'}, 'After', 'File');

    otherwise
        error('Unknown data_source: %s', string(data_source));
end
end

function T = make_lzc_file_table(csv_root, data_source)
switch string(data_source)
    case "full"
        p = fullfile(csv_root, 'PUBLIC_full_recording_LZc_by_file.csv');
        if ~exist(p, 'file')
            error(['Missing full-recording all-channel CSV:\n%s\n\n' ...
                'Run run_lzc_v7_full_recording_supplement.m first.'], p);
        end
        T = readtable(p, 'TextType', 'string');
        T.DataSource = repmat("full", height(T), 1);

    case {"window10s", "window2s"}
        p = fullfile(csv_root, 'PUBLIC_long_LZc_by_window.csv');
        if ~exist(p, 'file')
            error('V7 public all-channel window CSV not found: %s', p);
        end
        W = readtable(p, 'TextType', 'string');

        if string(data_source) == "window10s"
            window_sec = 10;
        else
            window_sec = 2;
        end

        W = W(asnum_local(W.WindowLengthSec) == window_sec, :);
        group_vars = {'Participant','Session','Eyes','Epoch','File'};
        T = groupsummary(W, group_vars, 'mean', ...
            {'LZc','LZcN','RawLZc','BinaryShuffleMeanRawLZc','PhaseRawLZcMean', ...
            'PropOnes','NTransitions','StringLength','N_Channels'});
        T.N_Windows = T.GroupCount;
        T.GroupCount = [];
        T.Properties.VariableNames = erase(T.Properties.VariableNames, 'mean_');
        T.WindowLengthSec = repmat(double(window_sec), height(T), 1);
        T.DataSource = repmat(string(data_source), height(T), 1);
        T = movevars(T, {'WindowLengthSec','DataSource','N_Windows'}, 'After', 'File');

    otherwise
        error('Unknown data_source: %s', string(data_source));
end
end

function build_per_file_channel_figures(T, data_source, metric, metric_color_limits, out_dir, template_labels, template_chanlocs)
file_keys = unique(T(:, {'Participant','Session','Eyes','Epoch','File'}), 'rows', 'stable');

for i = 1:height(file_keys)
    subset = T( ...
        T.Participant == file_keys.Participant(i) & ...
        T.Session == file_keys.Session(i) & ...
        T.Eyes == file_keys.Eyes(i) & ...
        T.Epoch == file_keys.Epoch(i) & ...
        T.File == file_keys.File(i), :);

    participant_dir = fullfile(out_dir, char(file_keys.Participant(i)));
    mkdir_if_needed(participant_dir);

    tag = make_clean_tag(file_keys.Participant(i), file_keys.Session(i), ...
        file_keys.Eyes(i), file_keys.Epoch(i), file_keys.File(i));

    plot_channel_heatmap(subset.(char(metric)), subset.Channel, ...
        [tag ' ' char(data_source) ' ' char(metric)], ...
        fullfile(participant_dir, [tag '_' char(data_source) '_' char(metric) '_channel_heatmap.png']), ...
        metric_color_limits);

    if exist('topoplot', 'file') == 2 && ~isempty(template_chanlocs)
        [aligned_values, aligned_chanlocs] = align_average_to_template( ...
            subset, template_labels, template_chanlocs, char(metric));

        plot_lzc_topoplot(aligned_values, aligned_chanlocs, ...
            [tag ' ' char(data_source) ' ' char(metric)], ...
            fullfile(participant_dir, [tag '_' char(data_source) '_' char(metric) '_topoplot.png']), ...
            metric_color_limits);
    end
end
end

function diff_table = make_eo_minus_ec_channel_table(T, metric)
keys = unique(T(:, {'Participant','Session','Epoch','Channel'}), 'rows', 'stable');
diff_table = table();
value_col = char(metric);
diff_col = [char(metric) '_EO_minus_EC'];

for i = 1:height(keys)
    eo = T(T.Participant == keys.Participant(i) & T.Session == keys.Session(i) & ...
        T.Epoch == keys.Epoch(i) & T.Channel == keys.Channel(i) & T.Eyes == "EO", :);
    ec = T(T.Participant == keys.Participant(i) & T.Session == keys.Session(i) & ...
        T.Epoch == keys.Epoch(i) & T.Channel == keys.Channel(i) & T.Eyes == "EC", :);

    if isempty(eo) || isempty(ec)
        continue;
    end

    row = table(keys.Participant(i), keys.Session(i), keys.Epoch(i), keys.Channel(i), ...
        double(mean(eo.(value_col), 'omitnan') - mean(ec.(value_col), 'omitnan')), ...
        'VariableNames', {'Participant','Session','Epoch','Channel', diff_col});
    diff_table = [diff_table; row];
end
end

function diff_table = make_eo_minus_ec_file_table(T, metric)
keys = unique(T(:, {'Participant','Session','Epoch'}), 'rows', 'stable');
diff_table = table();
value_col = char(metric);
diff_col = [char(metric) '_EO_minus_EC'];

for i = 1:height(keys)
    eo = T(T.Participant == keys.Participant(i) & T.Session == keys.Session(i) & ...
        T.Epoch == keys.Epoch(i) & T.Eyes == "EO", :);
    ec = T(T.Participant == keys.Participant(i) & T.Session == keys.Session(i) & ...
        T.Epoch == keys.Epoch(i) & T.Eyes == "EC", :);

    if isempty(eo) || isempty(ec)
        continue;
    end

    row = table(keys.Participant(i), keys.Session(i), keys.Epoch(i), ...
        double(mean(eo.(value_col), 'omitnan')), ...
        double(mean(ec.(value_col), 'omitnan')), ...
        double(mean(eo.(value_col), 'omitnan') - mean(ec.(value_col), 'omitnan')), ...
        'VariableNames', {'Participant','Session','Epoch', ...
        [char(metric) '_EO'], [char(metric) '_EC'], diff_col});
    diff_table = [diff_table; row];
end
end

function group_table = make_group_average_table(T, metric, session_order, eyes_order)
group_table = table();
metric_col = char(metric);
mean_col = ['Mean_' metric_col];
sd_col = ['SD_' metric_col];

for s = 1:numel(session_order)
    for e = 1:numel(eyes_order)
        subset = T(T.Session == string(session_order{s}) & T.Eyes == string(eyes_order{e}), :);
        if isempty(subset)
            continue;
        end

        channels = unique(subset.Channel, 'stable');
        for ch = 1:numel(channels)
            ch_subset = subset(subset.Channel == channels(ch), :);
            participant_means = groupsummary(ch_subset, 'Participant', 'mean', metric_col);
            participant_value_col = ['mean_' metric_col];
            row = table(string(session_order{s}), string(eyes_order{e}), channels(ch), ...
                double(mean(participant_means.(participant_value_col), 'omitnan')), ...
                double(std(participant_means.(participant_value_col), 'omitnan')), ...
                double(height(participant_means)), ...
                'VariableNames', {'Session','Eyes','Channel', mean_col, sd_col, 'N_Participants'});
            group_table = [group_table; row];
        end
    end
end
end

function group_table = make_file_group_summary(diff_table, metric, session_order)
group_table = table();
diff_col = [char(metric) '_EO_minus_EC'];

for s = 1:numel(session_order)
    subset_session = diff_table(diff_table.Session == string(session_order{s}), :);
    if isempty(subset_session)
        continue;
    end

    epochs = unique(subset_session.Epoch, 'stable');
    for e = 1:numel(epochs)
        subset = subset_session(subset_session.Epoch == epochs(e), :);
        values = subset.(diff_col);
        row = table(string(session_order{s}), epochs(e), ...
            double(mean(values, 'omitnan')), double(std(values, 'omitnan')), ...
            double(sum(isfinite(values))), ...
            'VariableNames', {'Session','Epoch','Mean_EO_minus_EC','SD_EO_minus_EC','N_Participants'});
        group_table = [group_table; row];
    end
end
end

function build_group_channel_figures(group_table, data_source, metric, metric_color_limits, group_fig_dir, ...
    template_labels, template_chanlocs, session_order, session_labels, eyes_order, eyes_labels)
mean_col = ['Mean_' char(metric)];

for s = 1:numel(session_order)
    for e = 1:numel(eyes_order)
        subset = group_table(group_table.Session == string(session_order{s}) & ...
            group_table.Eyes == string(eyes_order{e}), :);
        if isempty(subset)
            continue;
        end

        tag = sprintf('GROUP_%s_%s_%s_%s', session_order{s}, eyes_order{e}, char(data_source), char(metric));
        plot_channel_heatmap(subset.(mean_col), subset.Channel, ...
            [tag ' group average'], ...
            fullfile(group_fig_dir, [tag '_average_channel_heatmap.png']), metric_color_limits);

        if exist('topoplot', 'file') == 2 && ~isempty(template_chanlocs)
            [aligned_values, aligned_chanlocs] = align_average_to_template( ...
                subset, template_labels, template_chanlocs, mean_col);
            plot_lzc_topoplot(aligned_values, aligned_chanlocs, ...
                [tag ' group average'], ...
                fullfile(group_fig_dir, [tag '_average_topoplot.png']), metric_color_limits);
        end
    end
end

plot_group_2x5_channel_heatmaps(group_table, session_order, session_labels, eyes_order, eyes_labels, ...
    char(metric), fullfile(group_fig_dir, ...
    ['GROUP_ALL_SESSIONS_2x5_' char(data_source) '_' char(metric) '_channel_heatmaps.png']), ...
    metric_color_limits);

if exist('topoplot', 'file') == 2 && ~isempty(template_chanlocs)
    plot_group_2x5_topoplots(group_table, template_labels, template_chanlocs, ...
        session_order, session_labels, eyes_order, eyes_labels, ...
        char(metric), fullfile(group_fig_dir, ...
        ['GROUP_ALL_SESSIONS_2x5_' char(data_source) '_' char(metric) '_topoplots.png']), ...
        metric_color_limits);
end
end

function plot_file_group_summary(group_table, data_source, metric, out_fig)
if isempty(group_table)
    return;
end

labels = string(group_table.Session) + "_" + string(group_table.Epoch);
fig = figure('Color', 'w', 'Position', [100 100 1100 500]);
bar(group_table.Mean_EO_minus_EC);
hold on;
errorbar(1:height(group_table), group_table.Mean_EO_minus_EC, group_table.SD_EO_minus_EC, ...
    'k', 'LineStyle', 'none');
yline(0, '--k');
xticks(1:height(group_table));
xticklabels(labels);
xtickangle(45);
ylabel([char(metric) ' EO - EC']);
title(['All-channel ' char(metric) ' EO - EC summary (' char(data_source) ')'], 'Interpreter', 'none');
grid on;
saveas(fig, fullfile(out_fig, ['GROUP_EO_minus_EC_' char(data_source) '_' char(metric) '_summary.png']));
close(fig);
end

function [chanlocs, labels] = load_template_chanlocs(input_root)
chanlocs = [];
labels = strings(0, 1);
set_files = dir(fullfile(input_root, '**', 'PostICA_*.set'));

for i = 1:numel(set_files)
    try
        EEG = pop_loadset('filename', set_files(i).name, 'filepath', set_files(i).folder);
        if isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs)
            chanlocs = EEG.chanlocs;
            labels = get_channel_labels(EEG, EEG.nbchan);
            fprintf('Loaded template channel locations from:\n%s\n', fullfile(set_files(i).folder, set_files(i).name));
            return;
        end
    catch ME
        warning('Could not load template chanlocs from %s: %s', set_files(i).name, ME.message);
    end
end
warning('No template channel locations found. Topoplots will be skipped.');
end

function lims = percentile_limits(values, pct)
values = values(isfinite(values));
if isempty(values)
    lims = [0 1];
    return;
end
lims = prctile(values, pct);
if lims(1) == lims(2)
    pad = max(abs(lims(1)) * 0.05, 0.01);
    lims = [lims(1) - pad, lims(2) + pad];
end
end

function tag = make_clean_tag(participant, session_name, eyes, epoch, file)
tag = char(normalize_participant_id(participant) + "_" + string(session_name) + "_" + ...
    string(eyes) + "_" + string(epoch) + "_" + string(file));
tag = regexprep(tag, '[^\w-]', '_');
end

function x = asnum_local(x)
if isnumeric(x)
    return;
end
x = str2double(string(x));
end
