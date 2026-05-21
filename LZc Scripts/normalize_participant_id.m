function participant_id = normalize_participant_id(x)

if isnumeric(x)
    participant_id = string(sprintf('%03d', x));

elseif isstring(x)
    participant_id = x;

    if numel(participant_id) > 1
        participant_id = participant_id(1);
    end

    participant_char = char(participant_id);

    if ~isempty(participant_char) && all(isstrprop(participant_char, 'digit'))
        participant_id = string(sprintf('%03d', str2double(participant_char)));
    end

elseif iscell(x)
    participant_id = normalize_participant_id(x{1});

elseif ischar(x)
    participant_id = string(x);

    if ~isempty(x) && all(isstrprop(x, 'digit'))
        participant_id = string(sprintf('%03d', str2double(x)));
    end

else
    participant_id = string(x);
end

end