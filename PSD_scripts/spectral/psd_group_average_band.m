function group_band = psd_group_average_band(T)

T = normalize_psd_table_ids(T);

group_vars = {'Session','Eyes','Epoch','Channel','Band'};

group_band = groupsummary(T, group_vars, {'mean','std'}, ...
    {'AbsolutePower','Log10Power','RelativePower'});

end