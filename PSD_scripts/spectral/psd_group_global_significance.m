function stats_table = psd_group_global_significance(T, value_column, group_vars)
% Participant-level EO vs EC significance after averaging across channels.
% For band tables, use group_vars = {'Session','Epoch','Band'}.
% For aperiodic tables, use group_vars = {'Session','Epoch'}.

T = normalize_psd_table_ids(T);

required_vars = [group_vars, {'Participant','Eyes','Channel',value_column}];
for i = 1:numel(required_vars)
    if ~ismember(required_vars{i}, T.Properties.VariableNames)
        error('Missing required variable: %s', required_vars{i});
    end
end

keys_table = unique(T(:, group_vars), 'rows');
stats_table = table();

for k = 1:height(keys_table)

    subset = T;
    for g = 1:numel(group_vars)
        varname = group_vars{g};
        subset = subset(subset.(varname) == keys_table.(varname)(k), :);
    end

    if isempty(subset)
        continue;
    end

    % Average across channels first, separately for participant and eyes.
    P = groupsummary(subset, {'Participant','Eyes'}, 'mean', value_column);
    mean_col = ['mean_' value_column];

    eo = P(P.Eyes == "EO", {'Participant', mean_col});
    ec = P(P.Eyes == "EC", {'Participant', mean_col});

    if isempty(eo) || isempty(ec)
        continue;
    end

    joined = innerjoin(eo, ec, 'Keys', 'Participant');

    eo_col = [mean_col '_eo'];
    ec_col = [mean_col '_ec'];

    if ~ismember(eo_col, joined.Properties.VariableNames) || ~ismember(ec_col, joined.Properties.VariableNames)
        continue;
    end

    eo_vals = joined.(eo_col);
    ec_vals = joined.(ec_col);

    valid = isfinite(eo_vals) & isfinite(ec_vals);
    eo_vals = eo_vals(valid);
    ec_vals = ec_vals(valid);

    n = numel(eo_vals);
    if n < 3
        continue;
    end

    diff_vals = eo_vals - ec_vals;
    [~, p, ~, stats] = ttest(eo_vals, ec_vals);

    mean_diff = mean(diff_vals, 'omitnan');
    sd_diff = std(diff_vals, 'omitnan');

    row = keys_table(k,:);
    row.ValueColumn = string(value_column);
    row.N = n;
    row.Mean_EO_minus_EC = mean_diff;
    row.SD_EO_minus_EC = sd_diff;
    row.Tstat = stats.tstat;
    row.DF = stats.df;
    row.P = p;

    if sd_diff > 0
        row.Cohen_dz = mean_diff / sd_diff;
    else
        row.Cohen_dz = NaN;
    end

    stats_table = [stats_table; row];

end

if ~isempty(stats_table)
    stats_table.P_FDR = fdr_bh_psd(stats_table.P);
    stats_table.Significant_FDR_05 = stats_table.P_FDR < 0.05;
    stars = strings(height(stats_table),1);
    for i = 1:height(stats_table)
        stars(i) = psd_sig_star(stats_table.P_FDR(i));
    end
    stats_table.Significance = stars;
end

end
