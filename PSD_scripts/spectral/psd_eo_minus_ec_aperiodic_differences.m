function aperiodic_diff = psd_eo_minus_ec_aperiodic_differences(T)

T = normalize_psd_table_ids(T);

keys = {'Participant','Session','Epoch','Channel'};

eo = T(T.Eyes == "EO", :);
ec = T(T.Eyes == "EC", :);

eo.Eyes = [];
ec.Eyes = [];

joined = innerjoin(eo, ec, 'Keys', keys);

aperiodic_diff = joined(:, keys);

aperiodic_diff.AperiodicSlope_EO_minus_EC = joined.AperiodicSlope_eo - joined.AperiodicSlope_ec;
aperiodic_diff.AperiodicIntercept_EO_minus_EC = joined.AperiodicIntercept_eo - joined.AperiodicIntercept_ec;
aperiodic_diff.AperiodicRSquared_EO_minus_EC = joined.AperiodicRSquared_eo - joined.AperiodicRSquared_ec;

aperiodic_diff = normalize_psd_table_ids(aperiodic_diff);

end