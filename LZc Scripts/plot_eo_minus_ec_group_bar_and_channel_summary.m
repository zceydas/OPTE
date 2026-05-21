function plot_eo_minus_ec_group_bar_and_channel_summary(diff_results, plot_lzc_column, save_dir)

value_column = [char(plot_lzc_column) '_EO_minus_EC'];

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

target_labels = {
    'Baseline'
    'Dosing Epoch1'
    'Dosing Epoch2'
    'Dosing Epoch3'
    'Dosing Epoch4'
};

target_sessions = [
    "baseline"
    "dosing"
    "dosing"
    "dosing"
    "dosing"
];

target_epochs = [
    "Epoch0"
    "Epoch1"
    "Epoch2"
    "Epoch3"
    "Epoch4"
];

summary_rows = table();
channel_summary = table();

for i = 1:numel(target_labels)

    session_name = target_sessions(i);
    epoch_name = target_epochs(i);

    subset = diff_results( ...
        diff_results.Session == session_name & ...
        diff_results.Epoch == epoch_name, :);

    if isempty(subset)
        fprintf('No EO-EC data for %s %s. Skipping.\n', session_name, epoch_name);
        continue;
    end

    % Participant-level mean across all channels
    participant_means = groupsummary(subset, 'Participant', 'mean', value_column);

    group_mean = mean(participant_means.(['mean_' value_column]), 'omitnan');
    group_sd = std(participant_means.(['mean_' value_column]), 'omitnan');
    n_participants = height(participant_means);

    new_row = table( ...
        string(target_labels{i}), ...
        session_name, ...
        epoch_name, ...
        group_mean, ...
        group_sd, ...
        n_participants, ...
        'VariableNames', {'Condition','Session','Epoch','Mean_EO_minus_EC','SD_EO_minus_EC','N_Participants'} ...
    );

    summary_rows = [summary_rows; new_row];

    % Channel-level group mean across participants
    channel_means = groupsummary(subset, 'Channel', {'mean','std'}, value_column);

    condition_col = repmat(string(target_labels{i}), height(channel_means), 1);
    session_col = repmat(session_name, height(channel_means), 1);
    epoch_col = repmat(epoch_name, height(channel_means), 1);

    channel_table = table( ...
        condition_col, ...
        session_col, ...
        epoch_col, ...
        string(channel_means.Channel), ...
        channel_means.(['mean_' value_column]), ...
        channel_means.(['std_' value_column]), ...
        'VariableNames', {'Condition','Session','Epoch','Channel','Mean_EO_minus_EC','SD_EO_minus_EC'} ...
    );

    channel_summary = [channel_summary; channel_table];

end

summary_csv = fullfile(save_dir, ...
    ['EO_minus_EC_group_bar_summary_' char(plot_lzc_column) '.csv']);

channel_csv = fullfile(save_dir, ...
    ['EO_minus_EC_channel_summary_' char(plot_lzc_column) '.csv']);

writetable(summary_rows, summary_csv);
writetable(channel_summary, channel_csv);

%% Bar graph

fig1 = figure('Color', 'k', 'Position', [100 100 1000 600]);

bar(summary_rows.Mean_EO_minus_EC);
hold on;

errorbar( ...
    1:height(summary_rows), ...
    summary_rows.Mean_EO_minus_EC, ...
    summary_rows.SD_EO_minus_EC, ...
    'k', ...
    'LineStyle', 'none', ...
    'LineWidth', 1.5);

yline(0, '--k');

xticks(1:height(summary_rows));
xticklabels(summary_rows.Condition);
xtickangle(30);

ylabel('Mean EO - EC LZc across channels');
title(['Group EO - EC LZc Difference Summary using ' char(plot_lzc_column)], ...
    'Interpreter', 'none');

grid on;

bar_path = fullfile(save_dir, ...
    ['EO_minus_EC_group_bar_summary_' char(plot_lzc_column) '.png']);

saveas(fig1, bar_path);
close(fig1);

%% Channel × condition heatmap

channels = unique(channel_summary.Channel, 'stable');
conditions = string(target_labels);

heatmat = nan(numel(channels), numel(conditions));

for c = 1:numel(conditions)

    cond_subset = channel_summary(channel_summary.Condition == conditions(c), :);

    for ch = 1:numel(channels)
        idx = cond_subset.Channel == channels(ch);

        if any(idx)
            heatmat(ch, c) = cond_subset.Mean_EO_minus_EC(find(idx, 1));
        end
    end
end

fig2 = figure('Color', 'k', 'Position', [100 100 1000 1200]);

imagesc(heatmat);
colormap(parula);
colorbar;
caxis([-0.15 0.15]);

xticks(1:numel(conditions));
xticklabels(conditions);
xtickangle(30);

yticks(1:numel(channels));
yticklabels(channels);

xlabel('Condition');
ylabel('Channel');
title(['Channel-wise EO - EC LZc Difference using ' char(plot_lzc_column)], ...
    'Interpreter', 'none');

heatmap_path = fullfile(save_dir, ...
    ['EO_minus_EC_channel_by_condition_heatmap_' char(plot_lzc_column) '.png']);

saveas(fig2, heatmap_path);
close(fig2);

fprintf('\nSaved EO minus EC bar summary:\n%s\n', bar_path);
fprintf('Saved EO minus EC channel heatmap:\n%s\n', heatmap_path);
fprintf('Saved summary CSV:\n%s\n', summary_csv);
fprintf('Saved channel CSV:\n%s\n', channel_csv);

end