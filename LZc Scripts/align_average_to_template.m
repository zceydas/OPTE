function [aligned_values, aligned_chanlocs] = align_average_to_template(avg_table, template_labels, template_chanlocs, value_column)

aligned_values = nan(numel(template_labels), 1);

for i = 1:numel(template_labels)

    idx = avg_table.Channel == template_labels(i);

    if any(idx)
        aligned_values(i) = avg_table.(value_column)(find(idx, 1));
    end

end

valid_idx = ~isnan(aligned_values);

aligned_values = aligned_values(valid_idx);
aligned_chanlocs = template_chanlocs(valid_idx);

end