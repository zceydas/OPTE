function T = normalize_psd_table_ids(T)

if isempty(T)
    return;
end

if ismember('Participant', T.Properties.VariableNames)

    participant_clean = strings(height(T), 1);

    for i = 1:height(T)
        participant_clean(i) = normalize_participant_id(T.Participant(i));
    end

    T.Participant = participant_clean;
end

string_cols = {'Session','Eyes','Epoch','File','Channel','Band','Condition'};

for c = 1:numel(string_cols)
    col = string_cols{c};

    if ismember(col, T.Properties.VariableNames)
        T.(col) = string(T.(col));
    end
end

end