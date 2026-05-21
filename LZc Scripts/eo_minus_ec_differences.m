function diff_results = eo_minus_ec_differences(all_results, plot_lzc_column)

value_column = char(plot_lzc_column);

group_vars = {'Participant','Session','Eyes','Epoch','Channel'};

mean_table = groupsummary(all_results, group_vars, 'mean', value_column);

eo_table = mean_table(mean_table.Eyes == "EO", :);
ec_table = mean_table(mean_table.Eyes == "EC", :);

eo_table.Eyes = [];
ec_table.Eyes = [];

join_keys = {'Participant','Session','Epoch','Channel'};

joined = innerjoin(eo_table, ec_table, ...
    'Keys', join_keys);

left_col = ['mean_' value_column '_eo_table'];
right_col = ['mean_' value_column '_ec_table'];

if ~ismember(left_col, joined.Properties.VariableNames) || ...
        ~ismember(right_col, joined.Properties.VariableNames)

    possible = joined.Properties.VariableNames(contains(joined.Properties.VariableNames, ['mean_' value_column]));

    left_col = possible{1};
    right_col = possible{2};

end

diff_results = joined(:, join_keys);
diff_results.([value_column '_EO_minus_EC']) = joined.(left_col) - joined.(right_col);

end