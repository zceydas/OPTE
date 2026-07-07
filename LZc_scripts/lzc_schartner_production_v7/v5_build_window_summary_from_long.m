function window_summary = v5_build_window_summary_from_long(channel_long, lzc_long)
% Build window summary from long tables ONLY.

if isempty(lzc_long)
    window_summary = table();
    return;
end

keys = unique(lzc_long(:, {'Participant','Session','Eyes','Epoch','File','WindowSec','Window','WindowStartSec','WindowEndSec'}), 'rows');

window_summary = table();

for k = 1:height(keys)

    key = keys(k,:);

    ch_subset = channel_long( ...
        channel_long.Participant == key.Participant & ...
        channel_long.Session == key.Session & ...
        channel_long.Eyes == key.Eyes & ...
        channel_long.Epoch == key.Epoch & ...
        channel_long.File == key.File & ...
        channel_long.WindowSec == key.WindowSec & ...
        channel_long.Window == key.Window, :);

    lzc_subset = lzc_long( ...
        lzc_long.Participant == key.Participant & ...
        lzc_long.Session == key.Session & ...
        lzc_long.Eyes == key.Eyes & ...
        lzc_long.Epoch == key.Epoch & ...
        lzc_long.File == key.File & ...
        lzc_long.WindowSec == key.WindowSec & ...
        lzc_long.Window == key.Window, :);

    if height(lzc_subset) ~= 1
        error('Expected exactly one LZc row per window, found %d.', height(lzc_subset));
    end

    row = key;

    row.Mean_LZs = double(mean(asnum(ch_subset.LZs), 'omitnan'));
    row.SD_LZs = double(std(asnum(ch_subset.LZs), 'omitnan'));
    row.Mean_LZsN = double(mean(asnum(ch_subset.LZsN), 'omitnan'));
    row.SD_LZsN = double(std(asnum(ch_subset.LZsN), 'omitnan'));

    row.LZc = double(lzc_subset.LZc(1));
    row.LZcN = double(lzc_subset.LZcN(1));
    row.RawLZc = double(lzc_subset.RawLZc(1));
    row.BinaryShuffleMeanRawLZc = double(lzc_subset.BinaryShuffleMeanRawLZc(1));
    row.PhaseRawLZcMean = double(lzc_subset.PhaseRawLZcMean(1));
    row.N_Channels = double(lzc_subset.N_Channels(1));
    row.StringLength = double(lzc_subset.StringLength(1));

    window_summary = [window_summary; row];

end

end
