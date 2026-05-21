function plot_participant_dosing_2x4( ...
    all_results, participant_id, template_labels, template_chanlocs, ...
    eyes_order, eyes_labels, plot_lzc_column, lzc_color_limits, save_dir)

participant_id = normalize_participant_id(participant_id);
all_results = normalize_lzc_table_ids(all_results);

value_column = char(plot_lzc_column);

subset = all_results( ...
    all_results.Participant == participant_id & ...
    all_results.Session == "dosing", :);

if isempty(subset)
    return;
end

if isempty(template_chanlocs) || isempty(template_labels)
    warning('No template channel locations available. Skipping dosing 2x4 map.');
    return;
end

epochs = ["Epoch1","Epoch2","Epoch3","Epoch4"];

fig = figure('Color', 'w', 'Position', [50 50 1600 800]);
set(fig, 'InvertHardcopy', 'off');

tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for e = 1:2

    eyes_name = string(eyes_order{e});

    for ep = 1:4

        nexttile;

        epoch_name = epochs(ep);

        tile_data = subset( ...
            subset.Eyes == eyes_name & ...
            subset.Epoch == epoch_name, :);

        if isempty(tile_data)
            axis off;
            title(sprintf('%s\n%s\nNo data', eyes_labels{e}, epoch_name), ...
                'Interpreter', 'none', 'Color', 'k');
            continue;
        end

        avg_table = groupsummary(tile_data, 'Channel', 'mean', value_column);

        [aligned_values, aligned_chanlocs] = align_average_to_template( ...
            avg_table, template_labels, template_chanlocs, ['mean_' value_column]);

        if isempty(aligned_values) || isempty(aligned_chanlocs)
            axis off;
            title(sprintf('%s\n%s\nNo aligned channels', eyes_labels{e}, epoch_name), ...
                'Interpreter', 'none', 'Color', 'k');
            continue;
        end

        topoplot(aligned_values, aligned_chanlocs, ...
            'electrodes', 'off', ...
            'maplimits', lzc_color_limits, ...
            'plotrad', 0.5, ...
            'headrad', 0.5, ...
            'intrad', 0.5, ...
            'conv', 'off');

        title(sprintf('%s\n%s', eyes_labels{e}, epoch_name), ...
            'Interpreter', 'none', 'Color', 'k');

    end
end

cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb, char(plot_lzc_column));
cb.Color = 'k';

sgtitle(['Participant ' char(participant_id) ' Dosing LZc 2x4 Summary'], ...
    'Interpreter', 'none', 'Color', 'k');

save_path = fullfile(save_dir, ...
    ['Participant_' char(participant_id) '_dosing_2x4_' char(plot_lzc_column) '.png']);

saveas(fig, save_path);
close(fig);

end