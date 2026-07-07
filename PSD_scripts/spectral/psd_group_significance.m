function stats_table = psd_group_significance(T, value_column, group_vars)
% Group-level EO vs EC significance testing
T = normalize_psd_table_ids(T);
keys_table = unique(T(:, group_vars), 'rows');
stats_table = table();
for k = 1:height(keys_table)
    subset = T;
    for g = 1:numel(group_vars)
        subset = subset(subset.(group_vars{g}) == keys_table.(group_vars{g})(k), :);
    end
    eo = subset(subset.Eyes == "EO", :);
    ec = subset(subset.Eyes == "EC", :);
    if isempty(eo) || isempty(ec), continue; end
    eo.Eyes = []; ec.Eyes = [];
    joined = innerjoin(eo, ec, 'Keys', 'Participant');
    eo_col = [value_column '_eo'];
    ec_col = [value_column '_ec'];
    if ~ismember(eo_col, joined.Properties.VariableNames) || ~ismember(ec_col, joined.Properties.VariableNames)
        continue;
    end
    a = joined.(eo_col); b = joined.(ec_col);
    valid = isfinite(a) & isfinite(b);
    a = a(valid); b = b(valid);
    if numel(a) < 3, continue; end
    [~,p,~,s] = ttest(a,b);
    d = a-b;
    row = keys_table(k,:);
    row.N = numel(a);
    row.P = p;
    row.Tstat = s.tstat;
    row.DF = s.df;
    row.Cohen_dz = mean(d,'omitnan')/std(d,'omitnan');
    stats_table = [stats_table; row];
end
if ~isempty(stats_table)
    stats_table.P_FDR = fdr_bh_psd(stats_table.P);
end
end
