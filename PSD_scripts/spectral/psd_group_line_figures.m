function psd_group_line_figures(all_psd_results, save_dir)

if isempty(all_psd_results)
    warning('No PSD frequency table available for group line figures.');
    return;
end

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

T = normalize_psd_table_ids(all_psd_results);
conditions = psd_condition_order();
eyes_list = unique(T.Eyes, 'stable');
measures = {'AbsolutePower','Log10Power','RelativePower'};
ylabels = {'Absolute power','log10 power','Relative power'};

for m = 1:numel(measures)

    measure = measures{m};

    for e = 1:numel(eyes_list)

        eyes_name = string(eyes_list(e));
        fig_path = fullfile(save_dir, ['GROUP_PSD_' char(measure) '_' char(eyes_name) '_lineplot.png']);

        if exist(fig_path, 'file')
            fprintf('Group PSD line figure already exists. Skipping:\n%s\n', fig_path);
            continue;
        end

        fig = figure('Color','w','Position',[100 100 1100 650]);
        set(fig, 'InvertHardcopy','off');
        hold on;

        plotted_any = false;

        for c = 1:height(conditions)

            subset = T(T.Eyes == eyes_name & ...
                       T.Session == conditions.Session(c) & ...
                       T.Epoch == conditions.Epoch(c), :);

            if isempty(subset)
                continue;
            end

            P = groupsummary(subset, {'Participant','Frequency'}, 'mean', measure);
            freqs = unique(P.Frequency);
            mean_vals = nan(numel(freqs),1);
            sem_vals = nan(numel(freqs),1);

            colname = ['mean_' measure];

            for fi = 1:numel(freqs)
                vals = P.(colname)(P.Frequency == freqs(fi));
                mean_vals(fi) = mean(vals, 'omitnan');
                sem_vals(fi) = std(vals, 'omitnan') ./ sqrt(sum(~isnan(vals)));
            end

            plot(freqs, mean_vals, 'LineWidth', 1.5, 'DisplayName', char(conditions.Condition(c)));
            plotted_any = true;
        end

        if ~plotted_any
            close(fig);
            continue;
        end

        xlabel('Frequency (Hz)', 'Color','k');
        ylabel(ylabels{m}, 'Color','k');
        title(['Group PSD ' char(measure) ' - ' char(eyes_name)], 'Interpreter','none', 'Color','k');
        legend('Location','bestoutside', 'Interpreter','none');
        grid on;
        set(gca, 'Color','w', 'XColor','k', 'YColor','k');

        saveas(fig, fig_path);
        close(fig);
    end
end

end
