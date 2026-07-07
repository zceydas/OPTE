function Tpub = v5_public_file_summary_table(T)
% Public-facing file summary table.

if isempty(T)
    Tpub = table();
    return;
if ismember('WindowSec', Tpub.Properties.VariableNames)
    Tpub.WindowLengthSec = Tpub.WindowSec;
    Tpub = movevars(Tpub, 'WindowLengthSec', 'After', 'WindowSec');
end

end

keep = intersect({ ...
    'Participant','Session','Eyes','Epoch','File','WindowSec','N_Windows', ...
    'Mean_LZs','SD_LZs','Mean_LZsN','SD_LZsN', ...
    'Mean_LZc','SD_LZc','Mean_LZcN','SD_LZcN', ...
    'Mean_RawLZc','Mean_BinaryShuffleMeanRawLZc','Mean_PhaseRawLZc', ...
    'N_Channels','StringLength','N_Samples','SamplingRate', ...
    'N_BinaryShuffles','N_PhaseSurrogates','ReferenceStyle','PipelineVersion','DateProcessed'}, T.Properties.VariableNames, 'stable');

Tpub = T(:, keep);

if ismember('WindowSec', Tpub.Properties.VariableNames)
    Tpub.WindowLengthSec = Tpub.WindowSec;
    Tpub = movevars(Tpub, 'WindowLengthSec', 'After', 'WindowSec');
end

end
