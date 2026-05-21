function plot_group_2x5_topoplots(avg_results, template_labels, template_chanlocs, session_order, session_labels, eyes_order, eyes_labels, plot_lzc_column, save_path, lzc_color_limits)

value_column = sprintf('Mean_%s', plot_lzc_column);

fig = figure('Color', 'k', 'Position', [50 50 1800 850]);
set(fig, 'InvertHardcopy', 'off');

tiledlayout(numel(eyes_order), numel(session_order), ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for e = 1:numel(eyes_order)

    for s = 1:numel(session_order)

        nexttile;

        session_name = string(session_order{s});
        eyes_name = string(eyes_order{e});

        subset = avg_results( ...
            avg_results.Session == session_name & ...
            avg_results.Eyes == eyes_name, :);

        if isempty(subset)
            axis off;
            title(sprintf('%s\n%s\nNo data', session_labels{s}, eyes_labels{e}), ...
                'Interpreter', 'none', 'Color', 'k');
            continue;
        end

        [aligned_values, aligned_chanlocs] = align_average_to_template( ...
            subset, template_labels, template_chanlocs, value_column);

        if isempty(aligned_values)
            axis off;
            title(sprintf('%s\n%s\nNo aligned channels', session_labels{s}, eyes_labels{e}), ...
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

        title(sprintf('%s\n%s', session_labels{s}, eyes_labels{e}), ...
            'Interpreter', 'none', 'Color', 'k');

    end
end

cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb, ['Mean ' char(plot_lzc_column)]);

sgtitle(['Group Average LZc Topoplots Across Sessions using ' char(plot_lzc_column)], ...
    'Interpreter', 'none', 'Color', 'k');

saveas(fig, save_path);
close(fig);

end