function plot_participant_difference_maps( ...
    diff_results, participant_id, template_labels, template_chanlocs, ...
    session_order, session_labels, plot_lzc_column, diff_color_limits, save_dir)

participant_id = normalize_participant_id(participant_id);
diff_results = normalize_lzc_table_ids(diff_results);

value_column = [char(plot_lzc_column) '_EO_minus_EC'];

subset = diff_results(diff_results.Participant == participant_id, :);

if isempty(subset)
    return;
end

if isempty(template_chanlocs) || isempty(template_labels)
    warning('No template channel locations available. Skipping participant difference map.');
    return;
end

plot_sessions = [
    "baseline"
    "dosing"
    "dosing"
    "dosing"
    "dosing"
    "1week"
    "2week"
    "1month"
];

plot_epochs = [
    "Epoch0"
    "Epoch1"
    "Epoch2"
    "Epoch3"
    "Epoch4"
    "Epoch0"
    "Epoch0"
    "Epoch0"
];

plot_labels = {
    'Baseline'
    'Dosing Epoch1'
    'Dosing Epoch2'
    'Dosing Epoch3'
    'Dosing Epoch4'
    '1 Week'
    '2 Week'
    '1 Month'
};

fig = figure('Color', 'w', 'Position', [50 50 2200 650]);
set(fig, 'InvertHardcopy', 'off');

tiledlayout(1, numel(plot_sessions), ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for s = 1:numel(plot_sessions)

    nexttile;

    session_name = plot_sessions(s);
    epoch_name = plot_epochs(s);

    tile_data = subset( ...
        subset.Session == session_name & ...
        subset.Epoch == epoch_name, :);

    if isempty(tile_data)
        axis off;
        title(sprintf('%s\nNo data', plot_labels{s}), ...
            'Interpreter', 'none', 'Color', 'k');
        continue;
    end

    avg_table = groupsummary(tile_data, 'Channel', 'mean', value_column);

    [aligned_values, aligned_chanlocs] = align_average_to_template( ...
        avg_table, template_labels, template_chanlocs, ['mean_' value_column]);

    if isempty(aligned_values) || isempty(aligned_chanlocs)
        axis off;
        title(sprintf('%s\nNo aligned channels', plot_labels{s}), ...
            'Interpreter', 'none', 'Color', 'k');
        continue;
    end

    topoplot(aligned_values, aligned_chanlocs, ...
        'electrodes', 'off', ...
        'maplimits', diff_color_limits, ...
        'plotrad', 0.5, ...
        'headrad', 0.5, ...
        'intrad', 0.5, ...
        'conv', 'off');

    title(plot_labels{s}, 'Interpreter', 'none', 'Color', 'k');

end

cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb, 'EO - EC LZc');
cb.Color = 'k';

sgtitle(['Participant ' char(participant_id) ' EO minus EC LZc Difference'], ...
    'Interpreter', 'none', 'Color', 'k');

save_path = fullfile(save_dir, ...
    ['Participant_' char(participant_id) '_EO_minus_EC_' char(plot_lzc_column) '.png']);

saveas(fig, save_path);
close(fig);

end