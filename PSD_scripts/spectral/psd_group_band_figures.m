function psd_group_band_figures(all_band_results, save_dir, stats_abs, stats_log, stats_rel)

if nargin < 3, stats_abs = table(); end
if nargin < 4, stats_log = table(); end
if nargin < 5, stats_rel = table(); end

if isempty(all_band_results)
    warning('No band table available for group band figures.');
    return;
end

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

T = normalize_psd_table_ids(all_band_results);
conditions = psd_condition_order();
eyes_list = unique(T.Eyes, 'stable');
bands = unique(T.Band, 'stable');
measures = {'AbsolutePower','Log10Power','RelativePower'};
ylabels = {'Absolute power','log10 power','Relative power'};
stats_tables = {stats_abs, stats_log, stats_rel};

for m = 1:numel(measures)
    measure = measures{m};
    S = normalize_psd_table_ids(stats_tables{m});

    for e = 1:numel(eyes_list)
        eyes_name = string(eyes_list(e));
        fig_path = fullfile(save_dir, ['GROUP_band_' char(measure) '_' char(eyes_name) '_with_significance.png']);

        if exist(fig_path, 'file')
            fprintf('Group band figure already exists. Skipping:\n%s\n', fig_path);
            continue;
        end

        means = nan(height(conditions), numel(bands));
        sems = nan(height(conditions), numel(bands));
        pvals = nan(height(conditions), numel(bands));
        stars = strings(height(conditions), numel(bands));

        for c = 1:height(conditions)
            for b = 1:numel(bands)
                subset = T(T.Eyes == eyes_name & ...
                           T.Session == conditions.Session(c) & ...
                           T.Epoch == conditions.Epoch(c) & ...
                           T.Band == bands(b), :);

                if isempty(subset)
                    continue;
                end

                P = groupsummary(subset, 'Participant', 'mean', measure);
                colname = ['mean_' measure];
                vals = P.(colname);
                means(c,b) = mean(vals, 'omitnan');
                sems(c,b) = std(vals, 'omitnan') ./ sqrt(sum(~isnan(vals)));

                if ~isempty(S) && ismember('P_FDR', S.Properties.VariableNames)
                    idx = S.Session == conditions.Session(c) & ...
                          S.Epoch == conditions.Epoch(c) & ...
                          S.Band == bands(b);
                    if any(idx)
                        pvals(c,b) = S.P_FDR(find(idx,1));
                        stars(c,b) = psd_sig_star(pvals(c,b));
                    end
                end
            end
        end

        if all(isnan(means), 'all')
            continue;
        end

        fig = figure('Color','w','Position',[100 100 1200 650]);
        set(fig, 'InvertHardcopy','off');
        bh = bar(means, 'grouped');
        hold on;

        all_y = means(:) + sems(:);
        y_range = range(all_y(isfinite(all_y)));
        if isempty(y_range) || y_range == 0 || ~isfinite(y_range), y_range = 1; end
        y_offset = 0.04 * y_range;

        for b = 1:numel(bh)
            if isprop(bh(b), 'XEndPoints')
                x = bh(b).XEndPoints;
            else
                x = (1:height(conditions)) + (b - (numel(bh)+1)/2) * 0.12;
            end
            errorbar(x, means(:,b), sems(:,b), 'k', 'LineStyle','none', 'LineWidth',1);

            for c = 1:height(conditions)
                if strlength(stars(c,b)) > 0 && isfinite(means(c,b))
                    text(x(c), means(c,b) + sems(c,b) + y_offset, char(stars(c,b)), ...
                        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                        'Color','k', 'FontSize',14, 'FontWeight','bold');
                end
            end
        end

        xticks(1:height(conditions));
        xticklabels(conditions.Condition);
        xtickangle(30);
        ylabel(ylabels{m}, 'Color','k');
        title(['Group band power ' char(measure) ' - ' char(eyes_name) ' (* = EO vs EC, FDR < .05)'], ...
            'Interpreter','none', 'Color','k');
        legend(cellstr(bands), 'Location','bestoutside', 'Interpreter','none');
        grid on;
        set(gca, 'Color','w', 'XColor','k', 'YColor','k');

        saveas(fig, fig_path);
        close(fig);
    end
end

end
