%% spectral_analysis.m
clear; clc;
format long g

%% ---------------- PATH CONFIGURATION ----------------
% Set USE_HARDCODED_PATHS=true for reproducible batch runs.
% Leave it false to select folders interactively when the script runs.
USE_HARDCODED_PATHS = false;

% Optional hardcoded paths. Edit these for your machine if needed.
HARD_CODED_EEGLAB_PATH = '';    % e.g., '/path/to/eeglab2026.0.0'
HARD_CODED_INPUT_ROOT = '';     % folder containing PostICA_*.set files
HARD_CODED_OUTPUT_ROOT = '';    % folder where PSD_Results should be written

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

do_psd = true;
do_plot_psd = true;
do_spectrogram = true;
do_eo_minus_ec = true;
do_group_averages = true;
do_group_significance = true;
do_group_figures = true;

% Welch PSD settings
welch_window_sec = 4;
welch_overlap_fraction = 0.50;
nfft_min = 2048;

% Spectrogram settings
spectrogram_window_sec = 4;
spectrogram_step_sec = 1;
spectrogram_tapers = [3 5];   % [time-bandwidth product, number of tapers]

% Frequency settings
freq_min = 1;
freq_max = 45;
spectrogram_freqs = 1:0.5:45;

% Canonical frequency bands
bands = struct();
bands.delta = [1 4];
bands.theta = [4 8];
bands.alpha = [8 13];
bands.beta  = [13 30];
bands.gamma = [30 45];

band_names = fieldnames(bands);

% 1/f settings
aperiodic_fit_range = [2 40];
exclude_aperiodic_ranges = [
    17 18
    34 36
];

%% ---------------- SELECT INPUT FOLDER ----------------

if USE_HARDCODED_PATHS
    input_root = char(HARD_CODED_INPUT_ROOT);
    if isempty(input_root)
        error('HARD_CODED_INPUT_ROOT is empty. Set it or use interactive folder selection.');
    end
else
    input_root = uigetdir(pwd, 'Select folder containing preprocessed PostICA EEGLAB .set files');
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
    output_root = fullfile(input_root, 'PSD_Results');
end

if ~exist(output_root, 'dir')
    mkdir(output_root);
end

%% ---------------- FIND FILES ----------------

set_files = dir(fullfile(input_root, '**', 'PostICA_*.set'));

exclude_folders = string({set_files.folder})';
set_files = set_files(~contains(exclude_folders, 'PSD_Results'));
set_files = set_files(~contains(string({set_files.folder})', 'LZc_Results'));
set_files = set_files(~contains(string({set_files.folder})', 'Spectral_Results'));

if isempty(set_files)
    error('No PostICA_*.set files found under: %s', input_root);
end

fprintf('\nFound %d PostICA .set files.\n', numel(set_files));

%% ---------------- PARSER PREFLIGHT ----------------

parser_check = table();

for i = 1:numel(set_files)

    [~, base_name, ~] = fileparts(set_files(i).name);
    info = parse_psd_filename(base_name);

    parser_check = [parser_check; table( ...
        string(set_files(i).folder), ...
        string(set_files(i).name), ...
        string(info.participant), ...
        string(info.session), ...
        string(info.eyes), ...
        string(info.epoch), ...
        info.skip, ...
        'VariableNames', {'Folder','File','Participant','Session','Eyes','Epoch','Skip'} ...
    )];

end

writetable(parser_check, fullfile(output_root, 'PSD_parser_check.csv'));

fprintf('\nParser summary by participant:\n');
disp(groupsummary(parser_check, 'Participant'));

fprintf('\nParser summary by session/eyes/epoch:\n');
disp(groupsummary(parser_check, {'Session','Eyes','Epoch'}));

%% ---------------- MAIN LOOP ----------------

all_psd_results = table();
all_band_results = table();
all_aperiodic_results = table();
spectrogram_index = table();

for f = 1:numel(set_files)

    file = set_files(f).name;
    path = set_files(f).folder;
    set_path = fullfile(path, file);
    [~, base_name, ~] = fileparts(file);

    fprintf('\n====================================================\n');
    fprintf('Processing %d/%d:\n%s\n', f, numel(set_files), set_path);

    info = parse_psd_filename(base_name);

    if info.skip
        fprintf('Skipping file because parsing failed.\n');
        continue;
    end

    participant_dir = fullfile(output_root, char(info.participant));

    if ~exist(participant_dir, 'dir')
        mkdir(participant_dir);
    end

    file_tag = sprintf('%s_%s_%s_%s_%s', ...
        char(info.participant), char(info.session), char(info.eyes), char(info.epoch), base_name);

    psd_mat_path = fullfile(participant_dir, [file_tag '_PSD_outputs.mat']);
    psd_csv_path = fullfile(participant_dir, [file_tag '_PSD_by_channel_frequency.csv']);
    band_csv_path = fullfile(participant_dir, [file_tag '_band_power_by_channel.csv']);
    aperiodic_csv_path = fullfile(participant_dir, [file_tag '_aperiodic_1f_by_channel.csv']);

    psd_fig_path = fullfile(participant_dir, [file_tag '_PSD_lineplot.png']);
    spectrogram_mat_path = fullfile(participant_dir, [file_tag '_multitaper_spectrogram.mat']);
    spectrogram_fig_path = fullfile(participant_dir, [file_tag '_multitaper_spectrogram.png']);

    % Legacy fallback:
    % Older PSD runs may have been saved under a participant-specific PSD_Results
    % folder if the selected input folder was a single participant folder.
    legacy_output_root = fullfile(fileparts(path), 'PSD_Results');
    legacy_participant_dir = fullfile(legacy_output_root, char(info.participant));

    legacy_psd_mat_path = fullfile(legacy_participant_dir, [file_tag '_PSD_outputs.mat']);
    legacy_psd_csv_path = fullfile(legacy_participant_dir, [file_tag '_PSD_by_channel_frequency.csv']);
    legacy_band_csv_path = fullfile(legacy_participant_dir, [file_tag '_band_power_by_channel.csv']);
    legacy_aperiodic_csv_path = fullfile(legacy_participant_dir, [file_tag '_aperiodic_1f_by_channel.csv']);
    legacy_psd_fig_path = fullfile(legacy_participant_dir, [file_tag '_PSD_lineplot.png']);
    legacy_spectrogram_mat_path = fullfile(legacy_participant_dir, [file_tag '_multitaper_spectrogram.mat']);
    legacy_spectrogram_fig_path = fullfile(legacy_participant_dir, [file_tag '_multitaper_spectrogram.png']);

    if ~exist(psd_csv_path, 'file') && exist(legacy_psd_csv_path, 'file')
        fprintf('Using legacy PSD output location for this file:\n%s\n', legacy_participant_dir);
        psd_mat_path = legacy_psd_mat_path;
        psd_csv_path = legacy_psd_csv_path;
        band_csv_path = legacy_band_csv_path;
        aperiodic_csv_path = legacy_aperiodic_csv_path;
        psd_fig_path = legacy_psd_fig_path;
        spectrogram_mat_path = legacy_spectrogram_mat_path;
        spectrogram_fig_path = legacy_spectrogram_fig_path;
    end

    need_psd = ~exist(psd_mat_path, 'file') || ...
               ~exist(psd_csv_path, 'file') || ...
               ~exist(band_csv_path, 'file') || ...
               ~exist(aperiodic_csv_path, 'file');

    need_spectrogram = ~exist(spectrogram_mat_path, 'file') || ...
                       ~exist(spectrogram_fig_path, 'file');

    need_load = (do_psd && need_psd) || ...
            (do_spectrogram && need_spectrogram) || ...
            (do_plot_psd && ~exist(psd_fig_path, 'file'));

    if need_load
        EEG = pop_loadset('filename', file, 'filepath', path);
        X = double(EEG.data);
        [nChannels, ~] = size(X);
        chan_labels = get_channel_labels_psd(EEG, nChannels);
    end

    %% ---------------- PSD ----------------

    if do_psd && need_psd

        fprintf('Computing Welch PSD...\n');

        [psd_abs, freqs] = welch_psd_by_channel( ...
            X, EEG.srate, welch_window_sec, welch_overlap_fraction, ...
            nfft_min, freq_min, freq_max);

        psd_log = log10(psd_abs + eps);

        total_power = trapz(freqs, psd_abs, 2);
        psd_rel = psd_abs ./ total_power;

        psd_table = psd_to_long_table( ...
            psd_abs, psd_log, psd_rel, freqs, chan_labels, info, base_name);

        band_table = extract_band_power_table( ...
            psd_abs, psd_log, psd_rel, freqs, chan_labels, info, base_name, bands, band_names);

        aperiodic_table = aperiodic_1f_table( ...
            psd_abs, freqs, chan_labels, info, base_name, ...
            aperiodic_fit_range, exclude_aperiodic_ranges);

        writetable(psd_table, psd_csv_path);
        writetable(band_table, band_csv_path);
        writetable(aperiodic_table, aperiodic_csv_path);

        save(psd_mat_path, ...
            'psd_abs', 'psd_log', 'psd_rel', 'freqs', ...
            'chan_labels', 'band_table', 'aperiodic_table', ...
            'set_path', 'welch_window_sec', 'welch_overlap_fraction', ...
            'nfft_min', 'freq_min', 'freq_max', ...
            'bands', 'aperiodic_fit_range', 'exclude_aperiodic_ranges', ...
            '-v7.3');

    else
        if exist(psd_csv_path, 'file') && exist(band_csv_path, 'file') && ...
                exist(aperiodic_csv_path, 'file')

            fprintf('PSD outputs already exist. Loading CSVs/MAT.\n');
            psd_table = readtable(psd_csv_path, 'TextType', 'string');
            band_table = readtable(band_csv_path, 'TextType', 'string');
            aperiodic_table = readtable(aperiodic_csv_path, 'TextType', 'string');

            if exist(psd_mat_path, 'file')
                load(psd_mat_path, 'psd_abs', 'psd_log', 'psd_rel', 'freqs', 'chan_labels');
            end

        else
            if do_psd
                fprintf('Expected PSD outputs missing. Recomputing this file.\n');

                EEG = pop_loadset('filename', file, 'filepath', path);
                X = double(EEG.data);
                [nChannels, ~] = size(X);
                chan_labels = get_channel_labels_psd(EEG, nChannels);

                [psd_abs, freqs] = welch_psd_by_channel( ...
                    X, EEG.srate, welch_window_sec, welch_overlap_fraction, ...
                    nfft_min, freq_min, freq_max);

                psd_log = log10(psd_abs + eps);

                total_power = trapz(freqs, psd_abs, 2);
                psd_rel = psd_abs ./ total_power;

                psd_table = psd_to_long_table( ...
                    psd_abs, psd_log, psd_rel, freqs, chan_labels, info, base_name);

                band_table = extract_band_power_table( ...
                    psd_abs, psd_log, psd_rel, freqs, chan_labels, info, base_name, bands, band_names);

                aperiodic_table = aperiodic_1f_table( ...
                    psd_abs, freqs, chan_labels, info, base_name, ...
                    aperiodic_fit_range, exclude_aperiodic_ranges);

                writetable(psd_table, psd_csv_path);
                writetable(band_table, band_csv_path);
                writetable(aperiodic_table, aperiodic_csv_path);

                save(psd_mat_path, ...
                    'psd_abs', 'psd_log', 'psd_rel', 'freqs', ...
                    'chan_labels', 'band_table', 'aperiodic_table', ...
                    'set_path', 'welch_window_sec', 'welch_overlap_fraction', ...
                    'nfft_min', 'freq_min', 'freq_max', ...
                    'bands', 'aperiodic_fit_range', 'exclude_aperiodic_ranges', ...
                    '-v7.3');

            else
                warning(['PSD outputs are missing for this file, but do_psd is false. ' ...
                    'Skipping this file. Missing expected CSV: %s'], psd_csv_path);
                continue;
            end
        end
    end

    psd_table = normalize_psd_table_ids(psd_table);
    band_table = normalize_psd_table_ids(band_table);
    aperiodic_table = normalize_psd_table_ids(aperiodic_table);

    all_psd_results = [all_psd_results; psd_table];
    all_band_results = [all_band_results; band_table];
    all_aperiodic_results = [all_aperiodic_results; aperiodic_table];

    %% ---------------- PSD LINE PLOT ----------------

    if do_plot_psd && ~exist(psd_fig_path, 'file')

        if ~exist('psd_abs', 'var') || ~exist('freqs', 'var')
            load(psd_mat_path, 'psd_abs', 'psd_log', 'psd_rel', 'freqs', 'chan_labels');
        end

        plot_psd_lineplot( ...
            freqs, psd_abs, psd_log, psd_rel, ...
            [file_tag ' PSD'], psd_fig_path);

    else
        fprintf('PSD line plot already exists. Skipping.\n');
    end

    %% ---------------- MULTITAPER SPECTROGRAM ----------------

    if do_spectrogram && need_spectrogram

        fprintf('Computing multitaper spectrogram...\n');

        [spec_power, spec_times, spec_freqs] = multitaper_spectrogram_mean( ...
            X, EEG.srate, spectrogram_window_sec, spectrogram_step_sec, ...
            spectrogram_freqs, spectrogram_tapers);

        save(spectrogram_mat_path, ...
            'spec_power', 'spec_times', 'spec_freqs', ...
            'spectrogram_window_sec', 'spectrogram_step_sec', ...
            'spectrogram_freqs', 'spectrogram_tapers', ...
            'chan_labels', 'set_path', ...
            '-v7.3');

        plot_multitaper_spectrogram( ...
            spec_times, spec_freqs, spec_power, ...
            [file_tag ' Multitaper Spectrogram'], ...
            spectrogram_fig_path);

    else
        fprintf('Multitaper spectrogram already exists. Skipping.\n');
    end

    if exist(spectrogram_mat_path, 'file')
        spectrogram_index = [spectrogram_index; table( ...
            string(info.participant), string(info.session), string(info.eyes), string(info.epoch), ...
            string(base_name), string(spectrogram_mat_path), string(spectrogram_fig_path), ...
            'VariableNames', {'Participant','Session','Eyes','Epoch','File','SpectrogramMAT','SpectrogramPNG'} ...
        )];
    end

    clear psd_abs psd_log psd_rel freqs chan_labels spec_power spec_times spec_freqs psd_table band_table aperiodic_table

end

%% ---------------- SAVE COMBINED TABLES ----------------

all_psd_results = normalize_psd_table_ids(all_psd_results);
all_band_results = normalize_psd_table_ids(all_band_results);
all_aperiodic_results = normalize_psd_table_ids(all_aperiodic_results);
spectrogram_index = normalize_psd_table_ids(spectrogram_index);

combined_psd_csv = fullfile(output_root, 'ALL_PSD_by_channel_frequency.csv');
combined_band_csv = fullfile(output_root, 'ALL_PSD_band_power_by_channel.csv');
combined_aperiodic_csv = fullfile(output_root, 'ALL_PSD_aperiodic_1f_by_channel.csv');
spectrogram_index_csv = fullfile(output_root, 'ALL_multitaper_spectrogram_index.csv');

writetable(all_psd_results, combined_psd_csv);
writetable(all_band_results, combined_band_csv);
writetable(all_aperiodic_results, combined_aperiodic_csv);
writetable(spectrogram_index, spectrogram_index_csv);

save(fullfile(output_root, 'ALL_PSD_results.mat'), ...
    'all_psd_results', 'all_band_results', 'all_aperiodic_results', 'spectrogram_index', ...
    'bands', 'band_names', 'aperiodic_fit_range', 'exclude_aperiodic_ranges', ...
    '-v7.3');

fprintf('\nSaved combined PSD table:\n%s\n', combined_psd_csv);
fprintf('Saved combined PSD band table:\n%s\n', combined_band_csv);
fprintf('Saved combined 1/f table:\n%s\n', combined_aperiodic_csv);
fprintf('Saved spectrogram index:\n%s\n', spectrogram_index_csv);

%% ---------------- EO MINUS EC DIFFERENCES ----------------

band_diff = table();
aperiodic_diff = table();

if do_eo_minus_ec

    diff_dir = fullfile(output_root, 'EO_minus_EC_Differences');

    if ~exist(diff_dir, 'dir')
        mkdir(diff_dir);
    end

    band_diff_csv = fullfile(diff_dir, 'ALL_EO_minus_EC_PSD_band_power.csv');
    aperiodic_diff_csv = fullfile(diff_dir, 'ALL_EO_minus_EC_PSD_aperiodic_1f.csv');

    if exist(band_diff_csv, 'file') && exist(aperiodic_diff_csv, 'file')
        fprintf('\nEO minus EC PSD differences already exist. Loading.\n');
        band_diff = readtable(band_diff_csv, 'TextType', 'string');
        aperiodic_diff = readtable(aperiodic_diff_csv, 'TextType', 'string');
        band_diff = normalize_psd_table_ids(band_diff);
        aperiodic_diff = normalize_psd_table_ids(aperiodic_diff);
    else
        fprintf('\nComputing EO minus EC PSD differences.\n');

        band_diff = psd_eo_minus_ec_band_differences(all_band_results);
        aperiodic_diff = psd_eo_minus_ec_aperiodic_differences(all_aperiodic_results);

        writetable(band_diff, band_diff_csv);
        writetable(aperiodic_diff, aperiodic_diff_csv);
    end

end

%% ---------------- GROUP AVERAGES ----------------

if do_group_averages

    group_dir = fullfile(output_root, 'Group_Averages');

    if ~exist(group_dir, 'dir')
        mkdir(group_dir);
    end

    group_band = psd_group_average_band(all_band_results);
    group_aperiodic = psd_group_average_aperiodic(all_aperiodic_results);

    writetable(group_band, fullfile(group_dir, 'GROUP_PSD_band_power_by_channel.csv'));
    writetable(group_aperiodic, fullfile(group_dir, 'GROUP_PSD_aperiodic_1f_by_channel.csv'));

end


%% ---------------- GROUP SIGNIFICANCE ----------------

stats_abs = table();
stats_log = table();
stats_rel = table();
stats_slope = table();
stats_intercept = table();
stats_rsquared = table();

if do_group_significance

    stats_dir = fullfile(output_root, 'Group_Statistics');

    if ~exist(stats_dir, 'dir')
        mkdir(stats_dir);
    end

    fprintf('\nComputing group EO minus EC significance tests...\n');

    band_global_vars = {'Session','Epoch','Band'};

    stats_abs = psd_group_global_significance(all_band_results, 'AbsolutePower', band_global_vars);
    stats_log = psd_group_global_significance(all_band_results, 'Log10Power', band_global_vars);
    stats_rel = psd_group_global_significance(all_band_results, 'RelativePower', band_global_vars);

    writetable(stats_abs, fullfile(stats_dir, 'GROUP_global_stats_EO_minus_EC_absolute_power.csv'));
    writetable(stats_log, fullfile(stats_dir, 'GROUP_global_stats_EO_minus_EC_log10_power.csv'));
    writetable(stats_rel, fullfile(stats_dir, 'GROUP_global_stats_EO_minus_EC_relative_power.csv'));

    aperiodic_global_vars = {'Session','Epoch'};

    stats_slope = psd_group_global_significance(all_aperiodic_results, 'AperiodicSlope', aperiodic_global_vars);
    stats_intercept = psd_group_global_significance(all_aperiodic_results, 'AperiodicIntercept', aperiodic_global_vars);
    stats_rsquared = psd_group_global_significance(all_aperiodic_results, 'AperiodicRSquared', aperiodic_global_vars);

    writetable(stats_slope, fullfile(stats_dir, 'GROUP_global_stats_EO_minus_EC_aperiodic_slope.csv'));
    writetable(stats_intercept, fullfile(stats_dir, 'GROUP_global_stats_EO_minus_EC_aperiodic_intercept.csv'));
    writetable(stats_rsquared, fullfile(stats_dir, 'GROUP_global_stats_EO_minus_EC_aperiodic_rsquared.csv'));

    fprintf('Saved group significance tables to:\n%s\n', stats_dir);

end

%% ---------------- GROUP FIGURES ----------------

if do_group_figures

    group_fig_dir = fullfile(output_root, 'Group_Figures');

    if ~exist(group_fig_dir, 'dir')
        mkdir(group_fig_dir);
    end

    psd_group_line_figures(all_psd_results, group_fig_dir);
    psd_group_band_figures(all_band_results, group_fig_dir, stats_abs, stats_log, stats_rel);
    psd_group_aperiodic_figures(all_aperiodic_results, group_fig_dir, stats_slope, stats_intercept, stats_rsquared);

    if ~isempty(band_diff)
        psd_group_eo_minus_ec_band_figures(band_diff, group_fig_dir, stats_abs, stats_log, stats_rel);
    end

    if ~isempty(aperiodic_diff)
        psd_group_eo_minus_ec_aperiodic_figures(aperiodic_diff, group_fig_dir, stats_slope, stats_intercept, stats_rsquared);
    end

    if ~isempty(spectrogram_index)
        multitaper_group_spectrogram_figures(spectrogram_index, group_fig_dir);
    end

end

fprintf('\nPSD analysis complete.\n');
