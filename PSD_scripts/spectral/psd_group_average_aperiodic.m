function group_aperiodic = psd_group_average_aperiodic(T)

T = normalize_psd_table_ids(T);

group_vars = {'Session','Eyes','Epoch','Channel'};

group_aperiodic = groupsummary(T, group_vars, {'mean','std'}, ...
    {'AperiodicSlope','AperiodicIntercept','AperiodicRSquared'});

end