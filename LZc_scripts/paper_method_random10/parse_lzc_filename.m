function info = parse_lzc_filename(base_name)

name_lower = lower(base_name);

info = struct();
info.participant = "";
info.session = "";
info.eyes = "";
info.epoch = "EpochUnknown";
info.skip = false;

id_match = regexp(base_name, '(?<!\d)(\d{3,5})(?!\d)', 'tokens', 'once');

if isempty(id_match)
    info.skip = true;
    return;
else
    info.participant = string(id_match{1});
end

if contains(name_lower, 'baseline') || contains(name_lower, 'base') || ...
        contains(name_lower, 'ses-1') || contains(name_lower, 'session1') || ...
        contains(name_lower, 'tmseeg')
    info.session = "baseline";

elseif contains(name_lower, 'dosing') || contains(name_lower, 'dose') || ...
        contains(name_lower, 'ses-2') || contains(name_lower, 'session2')
    info.session = "dosing";

elseif contains(name_lower, '1week') || contains(name_lower, '1_week') || ...
        contains(name_lower, '1-week') || contains(name_lower, 'week1') || ...
        contains(name_lower, 'ses-3') || contains(name_lower, 'session3')
    info.session = "1week";

elseif contains(name_lower, '2week') || contains(name_lower, '2_week') || ...
        contains(name_lower, '2-week') || contains(name_lower, 'week2') || ...
        contains(name_lower, 'ses-4') || contains(name_lower, 'session4')
    info.session = "2week";

elseif contains(name_lower, '1month') || contains(name_lower, '1_month') || ...
        contains(name_lower, '1-month') || contains(name_lower, 'month1') || ...
        contains(name_lower, 'ses-5') || contains(name_lower, 'session5')
    info.session = "1month";

else
    info.skip = true;
    return;
end

if contains(name_lower, 'eyesopen') || contains(name_lower, 'eyes_open') || ...
        contains(name_lower, 'eyeopen') || contains(name_lower, '_eo') || ...
        contains(name_lower, '-eo') || contains(name_lower, ' eo')
    info.eyes = "EO";

elseif contains(name_lower, 'eyesclosed') || contains(name_lower, 'eyes_closed') || ...
        contains(name_lower, 'eyeclosed') || contains(name_lower, '_ec') || ...
        contains(name_lower, '-ec') || contains(name_lower, ' ec')
    info.eyes = "EC";

else
    info.skip = true;
    return;
end

epoch_match = regexp(base_name, 'Epoch[_-]?(\d+)', 'tokens', 'once', 'ignorecase');

if ~isempty(epoch_match)
    info.epoch = "Epoch" + string(epoch_match{1});
end

end