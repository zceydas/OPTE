function channel_summary = v5_build_channel_summary_from_long(channel_long)

if isempty(channel_long)
    channel_summary = table();
    return;
end

keys = unique(channel_long(:, {'Participant','Session','Eyes','Epoch','File','WindowSec','Channel','ChannelIndex'}), 'rows');

channel_summary = table();

for k = 1:height(keys)

    key = keys(k,:);

    subset = channel_long( ...
        channel_long.Participant == key.Participant & ...
        channel_long.Session == key.Session & ...
        channel_long.Eyes == key.Eyes & ...
        channel_long.Epoch == key.Epoch & ...
        channel_long.File == key.File & ...
        channel_long.WindowSec == key.WindowSec & ...
        channel_long.Channel == key.Channel, :);

    row = key;
    row.Mean_LZs = double(mean(asnum(subset.LZs), 'omitnan'));
    row.SD_LZs = double(std(asnum(subset.LZs), 'omitnan'));
    row.Mean_LZsN = double(mean(asnum(subset.LZsN), 'omitnan'));
    row.SD_LZsN = double(std(asnum(subset.LZsN), 'omitnan'));
    row.Mean_RawLZs = double(mean(asnum(subset.RawLZs), 'omitnan'));
    row.Mean_PhaseRawLZs = double(mean(asnum(subset.PhaseRawLZsMean), 'omitnan'));
    row.N_Windows = double(height(subset));

    channel_summary = [channel_summary; row];

end

end
