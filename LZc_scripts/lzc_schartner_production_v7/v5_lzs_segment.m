function result = v5_lzs_segment(x, seed, nBinaryShuffles, nPhaseSurrogates)

if nargin < 2; seed = []; end
if nargin < 3 || isempty(nBinaryShuffles); nBinaryShuffles = 1; end
if nargin < 4 || isempty(nPhaseSurrogates); nPhaseSurrogates = 10; end

[raw_lz, B, M, TH, binary_string] = v5_binary_string(x);

[~, shuf_mean, shuf_sd, shuf_values] = v5_binary_shuffle_norm(raw_lz, binary_string, seed, nBinaryShuffles);

if isfinite(shuf_mean) && shuf_mean ~= 0
    LZs = double(raw_lz / shuf_mean);
else
    LZs = NaN;
end

phase_raw = nan(nPhaseSurrogates, 1);

for s = 1:nPhaseSurrogates
    if isempty(seed)
        surr_seed = [];
    else
        surr_seed = seed + 10000 + s;
    end

    x_surr = v5_phase_shuffle_signal(x, surr_seed);
    [phase_raw(s), ~, ~, ~, ~] = v5_binary_string(x_surr);
end

phase_mean = double(mean(phase_raw, 'omitnan'));
phase_sd = double(std(phase_raw, 'omitnan'));

if isfinite(phase_mean) && phase_mean ~= 0
    LZsN = double(raw_lz / phase_mean);
else
    LZsN = NaN;
end

result = struct();
result.LZs = double(LZs);
result.LZsN = double(LZsN);
result.RawLZs = double(raw_lz);
result.BinaryShuffleMeanRawLZs = double(shuf_mean);
result.BinaryShuffleSDRawLZs = double(shuf_sd);
result.BinaryShuffleRawLZsValues = double(shuf_values(:));
result.PhaseRawLZsMean = double(phase_mean);
result.PhaseRawLZsSD = double(phase_sd);
result.PhaseRawLZsValues = double(phase_raw(:));
result.Threshold = double(TH(1));
result.PropOnes = double(mean(B(:)));
result.NTransitions = double(sum(diff(double(binary_string(:)')) ~= 0));
result.StringLength = double(numel(binary_string));
result.BinaryStringFirst1000 = string(binary_string(1:min(1000, numel(binary_string))));
result.B = B;
result.M = M;
result.TH = TH;

end
