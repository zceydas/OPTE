function plot_group_2x5_channel_heatmaps(avg_results, session_order, session_labels, eyes_order, eyes_labels, plot_lzc_column, save_path, lzc_color_limits)

value_column = sprintf('Mean_%s', plot_lzc_column);

fig = figure('Color', 'k', 'Position', [50 50 1800 700]);

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
                'Interpreter', 'none');
            continue;
        end

        imagesc(subset.(value_column)');
        caxis(lzc_color_limits);
        colormap(parula);

        yticks([]);
        xticks(1:height(subset));
        xticklabels(subset.Channel);
        xtickangle(90);

        title(sprintf('%s\n%s', session_labels{s}, eyes_labels{e}), ...
            'Interpreter', 'none');

        if e == numel(eyes_order)
            xlabel('Channel');
        end

    end
end

cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb, ['Mean ' char(plot_lzc_column)]);

sgtitle(['Group Average Channel-wise LZc Across Sessions using ' char(plot_lzc_column)], ...
    'Interpreter', 'none');

saveas(fig, save_path);
close(fig);

end