function psd_group_eo_minus_ec_aperiodic_figures(aperiodic_diff, save_dir, stats_slope, stats_intercept, stats_rsquared)

if nargin < 3, stats_slope = table(); end
if nargin < 4, stats_intercept = table(); end
if nargin < 5, stats_rsquared = table(); end

if isempty(aperiodic_diff)
    warning('No EO minus EC aperiodic table available for group figures.');
    return;
end

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

T = normalize_psd_table_ids(aperiodic_diff);

conditions = psd_condition_order();

metrics = {
    'AperiodicSlope_EO_minus_EC'
    'AperiodicIntercept_EO_minus_EC'
    'AperiodicRSquared_EO_minus_EC'
};

ylabels = {
    '1/f slope EO - EC'
    '1/f intercept EO - EC'
    '1/f R^2 EO - EC'
};

stats_tables = {stats_slope, stats_intercept, stats_rsquared};

fig_path = fullfile(save_dir, ...
    'GROUP_EO_minus_EC_aperiodic_1f_with_significance.png');

fig = figure('Color','w','Position',[100 100 1200 900]);
set(fig, 'InvertHardcopy','off');

tiledlayout(numel(metrics), 1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for m = 1:numel(metrics)

    metric = metrics{m};
    S = stats_tables{m};

    if ~isempty(S)
        S = normalize_psd_table_ids(S);
    end

    means = nan(height(conditions),1);
    sems = nan(height(conditions),1);
    stars = strings(height(conditions),1);

    for c = 1:height(conditions)

        subset = T( ...
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

            sig_idx = S.Session == conditions.Session(c) & ...
                      S.Epoch == conditions.Epoch(c);

            if any(sig_idx)
                stars(c) = psd_sig_star(S.P_FDR(find(sig_idx, 1)));
            end

        elseif ~isempty(S) && ismember('Significance', S.Properties.VariableNames)

            sig_idx = S.Session == conditions.Session(c) & ...
                      S.Epoch == conditions.Epoch(c);

            if any(sig_idx)
                stars(c) = string(S.Significance(find(sig_idx, 1)));
            end

        end

    end

    nexttile;

    bar(means);
    hold on;

    errorbar(1:height(conditions), means, sems, ...
        'k', 'LineStyle','none', 'LineWidth',1);

    yline(0, '--k');

    all_y = abs(means) + sems;
    y_range = range(all_y(isfinite(all_y)));

    if isempty(y_range) || y_range == 0 || ~isfinite(y_range)
        y_range = 1;
    end

    y_offset = 0.04 * y_range;

    for c = 1:height(conditions)

        if strlength(stars(c)) > 0 && isfinite(means(c))

            if means(c) >= 0
                y = means(c) + sems(c) + y_offset;
                valign = 'bottom';
            else
                y = means(c) - sems(c) - y_offset;
                valign = 'top';
            end

            text(c, y, char(stars(c)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment',valign, ...
                'Color','k', ...
                'FontSize',16, ...
                'FontWeight','bold');

        end
    end

    xticks(1:height(conditions));
    xticklabels(conditions.Condition);
    xtickangle(30);

    ylabel(ylabels{m}, 'Color','k');
    title(metric, 'Interpreter','none', 'Color','k');

    grid on;
    set(gca, 'Color','w', 'XColor','k', 'YColor','k');

end

sgtitle('Group EO - EC aperiodic 1/f differences (* = FDR corrected)', ...
    'Interpreter','none', 'Color','k');

saveas(fig, fig_path);
close(fig);

fprintf('Saved EO minus EC aperiodic figure with significance:\n%s\n', fig_path);

end