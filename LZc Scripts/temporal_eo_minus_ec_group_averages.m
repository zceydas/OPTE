function temporal_group_index = temporal_eo_minus_ec_group_averages( ...
    temporal_pair_index, session_order, session_labels, temporal_plot_lzc_column, ...
    diff_color_limits, temporal_group_dir)

temporal_group_index = table();

if isempty(temporal_pair_index)
    fprintf('No temporal pair differences available for group averaging.\n');
    return;
end

for s = 1:numel(session_order)

    session_name = string(session_order{s});
    session_label = session_labels{s};

    subset = temporal_pair_index(temporal_pair_index.Session == session_name, :);

    if isempty(subset)
        fprintf('No temporal EO minus EC differences for %s. Skipping.\n', session_name);
        continue;
    end

    group_mat = fullfile(temporal_group_dir, ...
        ['GROUP_' char(session_name) '_' char(temporal_plot_lzc_column) '_temporal_EO_minus_EC.mat']);

    group_png = fullfile(temporal_group_dir, ...
        ['GROUP_' char(session_name) '_' char(temporal_plot_lzc_column) '_temporal_EO_minus_EC_heatmap.png']);

    if exist(group_mat, 'file') && exist(group_png, 'file')
        fprintf('Group temporal EO minus EC already exists. Skipping:\n%s\n', group_mat);
    else
        loaded = cell(height(subset), 1);
        minCh = inf;
        minWin = inf;

        for i = 1:height(subset)
            temp = load(subset.Diff_MAT(i), 'temporal_diff', 'temporal_times_sec', 'chan_labels');
            loaded{i} = temp;

            minCh = min(minCh, size(temp.temporal_diff, 1));
            minWin = min(minWin, size(temp.temporal_diff, 2));
        end

        all_diffs = nan(minCh, minWin, height(subset));

        for i = 1:height(subset)
            all_diffs(:,:,i) = loaded{i}.temporal_diff(1:minCh, 1:minWin);
        end

        group_temporal_diff = mean(all_diffs, 3, 'omitnan');
        temporal_times_sec = loaded{1}.temporal_times_sec(1:minWin);
        chan_labels = loaded{1}.chan_labels(1:minCh);

        save(group_mat, ...
            'group_temporal_diff', ...
            'all_diffs', ...
            'temporal_times_sec', ...
            'chan_labels', ...
            'session_name', ...
            'session_label', ...
            'temporal_plot_lzc_column', ...
            'diff_color_limits', ...
            'subset', ...
            '-v7.3');

        temporal_lzc_heatmap( ...
            group_temporal_diff, temporal_times_sec, chan_labels, ...
            ['Group ' session_label ' Temporal EO - EC LZc'], ...
            group_png, ...
            diff_color_limits);

        fprintf('Saved group temporal EO minus EC heatmap:\n%s\n', group_png);
    end

    new_row = table( ...
        session_name, ...
        string(group_mat), ...
        string(group_png), ...
        height(subset), ...
        'VariableNames', {'Session','Group_MAT','Group_PNG','N_Pairs'} ...
    );

    temporal_group_index = [temporal_group_index; new_row];

end

end