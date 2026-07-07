function Tpub = v5_public_lzc_table(T)
% Public-facing long-format all-channel LZc per window table.

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
    'LZc','LZcN','RawLZc', ...
    'BinaryShuffleMeanRawLZc','BinaryShuffleSDRawLZc', ...
    'PhaseRawLZcMean','PhaseRawLZcSD', ...
    'PropOnes','NTransitions','StringLength','N_Channels'});

if ismember('WindowSec', Tpub.Properties.VariableNames)
    Tpub.WindowLengthSec = Tpub.WindowSec;
    Tpub = movevars(Tpub, 'WindowLengthSec', 'After', 'WindowSec');
end

end
