function Tpub = v5_public_channel_lzs_table(T)
% Public-facing long-format per-channel/window table.

if isempty(T)
    Tpub = table();
    return;
if ismember('WindowSec', Tpub.Properties.VariableNames)
    Tpub.WindowLengthSec = Tpub.WindowSec;
    Tpub = movevars(Tpub, 'WindowLengthSec', 'After', 'WindowSec');
end

end

Tpub = T(:, { ...
    'Participant','Session','Eyes','Epoch','File', ...
    'WindowSec','Window','WindowStartSec','WindowEndSec', ...
    'Channel','ChannelIndex', ...
    'LZs','LZsN','RawLZs', ...
    'BinaryShuffleMeanRawLZs','BinaryShuffleSDRawLZs', ...
    'PhaseRawLZsMean','PhaseRawLZsSD', ...
    'Threshold','PropOnes','NTransitions','StringLength'});

if ismember('WindowSec', Tpub.Properties.VariableNames)
    Tpub.WindowLengthSec = Tpub.WindowSec;
    Tpub = movevars(Tpub, 'WindowLengthSec', 'After', 'WindowSec');
end

end
