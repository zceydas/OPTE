function plot_lzc_topoplot(values, chanlocs, plot_title, save_path, lzc_color_limits)

valid_idx = ~isnan(values);

values = values(valid_idx);
chanlocs = chanlocs(valid_idx);

fig = figure('Color', 'k', 'Position', [100 100 700 600]);
set(fig, 'InvertHardcopy', 'off');

if isempty(values)
    axis off;
    title([plot_title ' - No valid values'], 'Interpreter', 'none', 'Color', 'k');
    saveas(fig, save_path);
    close(fig);
    return;
end

topoplot(values, chanlocs, ...
    'electrodes', 'off', ...
    'maplimits', lzc_color_limits, ...
    'plotrad', 0.5, ...
    'headrad', 0.5, ...
    'intrad', 0.5, ...
    'conv', 'off');

colorbar;
title(plot_title, 'Interpreter', 'none', 'Color', 'k');
set(gca, 'Color', 'k');

saveas(fig, save_path);
close(fig);

end