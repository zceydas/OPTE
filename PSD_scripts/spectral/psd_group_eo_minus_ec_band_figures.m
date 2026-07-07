function psd_group_eo_minus_ec_band_figures(band_diff, save_dir, stats_abs, stats_log, stats_rel)

if nargin < 3, stats_abs = table(); end
if nargin < 4, stats_log = table(); end
if nargin < 5, stats_rel = table(); end

if isempty(band_diff)
    warning('No EO minus EC band table available for group figures.');
    return;
end

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

T = normalize_psd_table_ids(band_diff);

conditions = psd_condition_order();
bands = unique(T.Band, 'stable');

measures = {
    'AbsolutePower_EO_minus_EC'
    'Log10Power_EO_minus_EC'
    'RelativePower_EO_minus_EC'
};

ylabels = {
    'Absolute power EO - EC'
    'log10 power EO - EC'
    'Relative power EO - EC'
};

stats_tables = {stats_abs, stats_log, stats_rel};

for m = 1:numel(measures)

    measure = measures{m};
    S = stats_tables{m};

    if ~isempty(S)
        S = normalize_psd_table_ids(S);
    end

    fig_path = fullfile(save_dir, ...
        ['GROUP_EO_minus_EC_band_' char(measure) '_with_significance.png']);

    means = nan(height(conditions), numel(bands));
    sems = nan(height(conditions), numel(bands));
    stars = strings(height(conditions), numel(bands));

    for c = 1:height(conditions)

        for b = 1:numel(bands)

            subset = T( ...
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

                sig_idx = S.Session == conditions.Session(c) & ...
                          S.Epoch == conditions.Epoch(c) & ...
                          S.Band == bands(b);

                if any(sig_idx)
                    stars(c,b) = psd_sig_star(S.P_FDR(find(sig_idx, 1)));
                end

            elseif ~isempty(S) && ismember('Significance', S.Properties.VariableNames)

                sig_idx = S.Session == conditions.Session(c) & ...
                          S.Epoch == conditions.Epoch(c) & ...
                          S.Band == bands(b);

                if any(sig_idx)
                    stars(c,b) = string(S.Significance(find(sig_idx, 1)));
                end

            end

        end
    end

    if all(isnan(means), 'all')
        continue;
    end

    fig = figure('Color','w','Position',[100 100 1300 700]);
    set(fig, 'InvertHardcopy','off');

    bh = bar(means, 'grouped');
    hold on;

    yline(0, '--k');

    all_y = abs(means(:)) + sems(:);
    y_range = range(all_y(isfinite(all_y)));

    if isempty(y_range) || y_range == 0 || ~isfinite(y_range)
        y_range = 1;
    end

    y_offset = 0.04 * y_range;

    for b = 1:numel(bh)

        if isprop(bh(b), 'XEndPoints')
            x = bh(b).XEndPoints;
        else
            x = (1:height(conditions)) + ...
                (b - (numel(bh)+1)/2) * 0.12;
        end

        errorbar(x, means(:,b), sems(:,b), ...
            'k', 'LineStyle','none', 'LineWidth',1);

        for c = 1:height(conditions)

            if strlength(stars(c,b)) > 0 && isfinite(means(c,b))

                if means(c,b) >= 0
                    y = means(c,b) + sems(c,b) + y_offset;
                    valign = 'bottom';
                else
                    y = means(c,b) - sems(c,b) - y_offset;
                    valign = 'top';
                end

                text(x(c), y, char(stars(c,b)), ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment',valign, ...
                    'Color','k', ...
                    'FontSize',16, ...
                    'FontWeight','bold');

            end
        end
    end

    xticks(1:height(conditions));
    xticklabels(conditions.Condition);
    xtickangle(30);

    ylabel(ylabels{m}, 'Color','k');
    title(['Group EO - EC band differences: ' char(measure) ' (* = FDR corrected)'], ...
        'Interpreter','none', 'Color','k');

    legend(cellstr(bands), ...
        'Location','bestoutside', ...
        'Interpreter','none');

    grid on;
    set(gca, 'Color','w', 'XColor','k', 'YColor','k');

    saveas(fig, fig_path);
    close(fig);

    fprintf('Saved EO minus EC band figure with significance:\n%s\n', fig_path);

end

end