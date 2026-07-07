function band_diff = psd_eo_minus_ec_band_differences(T)

T = normalize_psd_table_ids(T);

keys = {'Participant','Session','Epoch','Channel','Band'};

eo = T(T.Eyes == "EO", :);
ec = T(T.Eyes == "EC", :);

eo.Eyes = [];
ec.Eyes = [];

joined = innerjoin(eo, ec, 'Keys', keys);

band_diff = joined(:, keys);

band_diff.AbsolutePower_EO_minus_EC = joined.AbsolutePower_eo - joined.AbsolutePower_ec;
band_diff.Log10Power_EO_minus_EC = joined.Log10Power_eo - joined.Log10Power_ec;
band_diff.RelativePower_EO_minus_EC = joined.RelativePower_eo - joined.RelativePower_ec;

band_diff = normalize_psd_table_ids(band_diff);

end