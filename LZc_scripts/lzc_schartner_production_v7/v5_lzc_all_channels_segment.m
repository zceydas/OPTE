function result = v5_lzc_all_channels_segment(X, seed, nBinaryShuffles, nPhaseSurrogates)

if nargin < 2; seed = []; end
if nargin < 3 || isempty(nBinaryShuffles); nBinaryShuffles = 1; end
if nargin < 4 || isempty(nPhaseSurrogates); nPhaseSurrogates = 10; end

[raw_lz, B, M, TH, binary_string] = v5_binary_string(X);

[~, shuf_mean, shuf_sd, shuf_values] = v5_binary_shuffle_norm(raw_lz, binary_string, seed, nBinaryShuffles);

if isfinite(shuf_mean) && shuf_mean ~= 0
    LZc = double(raw_lz / shuf_mean);
else
    LZc = NaN;
end

phase_raw = nan(nPhaseSurrogates, 1);

for s = 1:nPhaseSurrogates
    if isempty(seed)
        surr_seed = [];
    else
        surr_seed = seed + 10000 + s;
    end

    X_surr = v5_phase_shuffle_multichannel(X, surr_seed);
    [phase_raw(s), ~, ~, ~, ~] = v5_binary_string(X_surr);
end

phase_mean = double(mean(phase_raw, 'omitnan'));
phase_sd = double(std(phase_raw, 'omitnan'));

if isfinite(phase_mean) && phase_mean ~= 0
    LZcN = double(raw_lz / phase_mean);
else
    LZcN = NaN;
end

result = struct();
result.LZc = double(LZc);
result.LZcN = double(LZcN);
result.RawLZc = double(raw_lz);
result.BinaryShuffleMeanRawLZc = double(shuf_mean);
result.BinaryShuffleSDRawLZc = double(shuf_sd);
result.BinaryShuffleRawLZcValues = double(shuf_values(:));
result.PhaseRawLZcMean = double(phase_mean);
result.PhaseRawLZcSD = double(phase_sd);
result.PhaseRawLZcValues = double(phase_raw(:));
result.PropOnes = double(mean(B(:)));
result.NTransitions = double(sum(diff(double(binary_string(:)')) ~= 0));
result.StringLength = double(numel(binary_string));
result.BinaryStringFirst1000 = string(binary_string(1:min(1000, numel(binary_string))));
result.B = B;
result.M = M;
result.TH = TH;

end
