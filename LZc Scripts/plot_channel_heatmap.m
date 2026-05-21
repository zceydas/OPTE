function plot_channel_heatmap(values, chan_labels, plot_title, save_path, lzc_color_limits)

fig = figure('Color', 'k', 'Position', [100 100 1300 450]);

imagesc(values');
colormap(parula);
colorbar;
caxis(lzc_color_limits);

yticks([]);
xticks(1:numel(values));
xticklabels(chan_labels);
xtickangle(90);

title(plot_title, 'Interpreter', 'none');
xlabel('Channel');

saveas(fig, save_path);
close(fig);

end