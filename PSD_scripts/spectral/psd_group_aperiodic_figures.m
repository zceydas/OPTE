function psd_group_aperiodic_figures(all_aperiodic_results, save_dir, stats_slope, stats_intercept, stats_rsquared)

if nargin < 3, stats_slope = table(); end
if nargin < 4, stats_intercept = table(); end
if nargin < 5, stats_rsquared = table(); end

if isempty(all_aperiodic_results)
    warning('No aperiodic table available for group 1/f figures.');
    return;
end

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

T = normalize_psd_table_ids(all_aperiodic_results);
conditions = psd_condition_order();
eyes_list = unique(T.Eyes, 'stable');
metrics = {'AperiodicSlope','AperiodicIntercept','AperiodicRSquared'};
ylabels = {'1/f slope','1/f intercept','1/f R^2'};
stats_tables = {stats_slope, stats_intercept, stats_rsquared};

for e = 1:numel(eyes_list)
    eyes_name = string(eyes_list(e));
    fig_path = fullfile(save_dir, ['GROUP_aperiodic_1f_' char(eyes_name) '_with_significance.png']);

    if exist(fig_path, 'file')
        fprintf('Group aperiodic figure already exists. Skipping:\n%s\n', fig_path);
        continue;
    end

    fig = figure('Color','w','Position',[100 100 1200 850]);
    set(fig, 'InvertHardcopy','off');
    tiledlayout(numel(metrics), 1, 'TileSpacing','compact', 'Padding','compact');

    for m = 1:numel(metrics)
        metric = metrics{m};
        S = normalize_psd_table_ids(stats_tables{m});
        means = nan(height(conditions),1);
        sems = nan(height(conditions),1);
        stars = strings(height(conditions),1);

        for c = 1:height(conditions)
            subset = T(T.Eyes == eyes_name & ...
                       T.Session == conditions.Session(c) & ...
                       T.Epoch == conditions.Epoch(c), :);

            if isempty(subset)
                continue;
            end

            P = groupsummary(subset, 'Participant', 'mean', metric);
            colname = ['mean_' metric];
            vals = P.(colname);
            means(c) = mean(vals, 'omitnan');
            sems(c) = std(vals, 'omitnan') ./ sqrt(sum(~isnan(vals)));

            if ~isempty(S) && ismember('P_FDR', S.Properties.VariableNames)
                idx = S.Session == conditions.Session(c) & S.Epoch == conditions.Epoch(c);
                if any(idx)
                    stars(c) = psd_sig_star(S.P_FDR(find(idx,1)));
                end
            end
        end

        nexttile;
        bar(means);
        hold on;
        errorbar(1:height(conditions), means, sems, 'k', 'LineStyle','none', 'LineWidth',1);
        yline(0, '--k');

        all_y = means + sems;
        y_range = range(all_y(isfinite(all_y)));
        if isempty(y_range) || y_range == 0 || ~isfinite(y_range), y_range = 1; end
        y_offset = 0.04 * y_range;

        for c = 1:height(conditions)
            if strlength(stars(c)) > 0 && isfinite(means(c))
                text(c, means(c) + sems(c) + y_offset, char(stars(c)), ...
                    'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                    'Color','k', 'FontSize',14, 'FontWeight','bold');
            end
        end

        xticks(1:height(conditions));
        xticklabels(conditions.Condition);
        xtickangle(30);
        ylabel(ylabels{m}, 'Color','k');
        title([char(metric) ' - ' char(eyes_name)], 'Interpreter','none', 'Color','k');
        grid on;
        set(gca, 'Color','w', 'XColor','k', 'YColor','k');
    end

    sgtitle(['Group aperiodic 1/f metrics - ' char(eyes_name) ' (* = EO vs EC, FDR < .05)'], ...
        'Interpreter','none', 'Color','k');
    saveas(fig, fig_path);
    close(fig);
end

end
