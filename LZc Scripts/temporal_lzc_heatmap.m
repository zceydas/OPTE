function temporal_lzc_heatmap(temporal_LZc, temporal_times_sec, chan_labels, plot_title, save_path, lzc_color_limits)

fig = figure('Color', 'w', 'Position', [100 100 1400 700]);

imagesc(temporal_times_sec, 1:numel(chan_labels), temporal_LZc);
set(gca, 'YDir', 'normal');

colormap(parula);
colorbar;
caxis(lzc_color_limits);

yticks(1:numel(chan_labels));
yticklabels(chan_labels);

xlabel('Time (seconds)');
ylabel('Channel');
title(plot_title, 'Interpreter', 'none');

saveas(fig, save_path);
close(fig);

end