function plot_multitaper_spectrogram(spec_times, spec_freqs, spec_power, plot_title, save_path)

fig = figure('Color','w','Position',[100 100 1200 600]);
set(fig, 'InvertHardcopy', 'off');

imagesc(spec_times, spec_freqs, log10(spec_power + eps));
axis xy;

xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(plot_title, 'Interpreter','none','Color','k');

cb = colorbar;
ylabel(cb, 'log10 power');
cb.Color = 'k';

set(gca, 'Color','w', 'XColor','k', 'YColor','k');

saveas(fig, save_path);
close(fig);

end