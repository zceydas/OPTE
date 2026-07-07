function file_summary = v5_build_file_summary_from_window_summary(window_summary, nSamples, srate, nBinaryShuffles, nPhaseSurrogates)

if isempty(window_summary)
    file_summary = table();
    return;
end

base = window_summary(1, {'Participant','Session','Eyes','Epoch','File','WindowSec'});

file_summary = base;
file_summary.N_Windows = double(height(window_summary));

file_summary.Mean_LZs = double(mean(asnum(window_summary.Mean_LZs), 'omitnan'));
file_summary.SD_LZs = double(std(asnum(window_summary.Mean_LZs), 'omitnan'));
file_summary.Mean_LZsN = double(mean(asnum(window_summary.Mean_LZsN), 'omitnan'));
file_summary.SD_LZsN = double(std(asnum(window_summary.Mean_LZsN), 'omitnan'));

file_summary.Mean_LZc = double(mean(asnum(window_summary.LZc), 'omitnan'));
file_summary.SD_LZc = double(std(asnum(window_summary.LZc), 'omitnan'));
file_summary.Mean_LZcN = double(mean(asnum(window_summary.LZcN), 'omitnan'));
file_summary.SD_LZcN = double(std(asnum(window_summary.LZcN), 'omitnan'));

file_summary.Mean_RawLZc = double(mean(asnum(window_summary.RawLZc), 'omitnan'));
file_summary.Mean_BinaryShuffleMeanRawLZc = double(mean(asnum(window_summary.BinaryShuffleMeanRawLZc), 'omitnan'));
file_summary.Mean_PhaseRawLZc = double(mean(asnum(window_summary.PhaseRawLZcMean), 'omitnan'));

file_summary.N_Channels = double(window_summary.N_Channels(1));
file_summary.StringLength = double(round(mean(asnum(window_summary.StringLength), 'omitnan')));

file_summary.N_Samples = double(nSamples);
file_summary.SamplingRate = double(srate);
file_summary.N_BinaryShuffles = double(nBinaryShuffles);
file_summary.N_PhaseSurrogates = double(nPhaseSurrogates);

end
