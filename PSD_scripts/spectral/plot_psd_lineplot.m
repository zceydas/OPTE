function plot_psd_lineplot(freqs, psd_abs, psd_log, psd_rel, plot_title, save_path)

mean_abs = mean(psd_abs, 1, 'omitnan');
mean_log = mean(psd_log, 1, 'omitnan');
mean_rel = mean(psd_rel, 1, 'omitnan');

fig = figure('Color','w','Position',[100 100 1200 700]);
set(fig, 'InvertHardcopy', 'off');

tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(freqs, mean_abs, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Absolute power');
title('Absolute PSD');
grid on;
set(gca, 'Color','w', 'XColor','k', 'YColor','k');

nexttile;
plot(freqs, mean_log, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('log10 power');
title('Log10 PSD');
grid on;
set(gca, 'Color','w', 'XColor','k', 'YColor','k');

nexttile;
plot(freqs, mean_rel, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Relative power');
title('Relative PSD');
grid on;
set(gca, 'Color','w', 'XColor','k', 'YColor','k');

sgtitle(plot_title, 'Interpreter','none','Color','k');

saveas(fig, save_path);
close(fig);

end