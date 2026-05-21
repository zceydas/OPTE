function [temporal_LZc, temporal_times_sec] = temporal_lzc_by_channel( ...
    X, srate, window_sec, step_sec, seed, nShuffles)

[nChannels, nSamples] = size(X);

window_samples = round(window_sec * srate);
step_samples = round(step_sec * srate);

start_samples = 1:step_samples:(nSamples - window_samples + 1);
nWindows = numel(start_samples);

temporal_LZc = nan(nChannels, nWindows);
temporal_times_sec = nan(1, nWindows);

for w = 1:nWindows

    idx1 = start_samples(w);
    idx2 = idx1 + window_samples - 1;

    temporal_times_sec(w) = ((idx1 + idx2) / 2) / srate;

    fprintf('    Temporal window %d/%d\n', w, nWindows);

    for ch = 1:nChannels

        X_seg = X(ch, idx1:idx2);

        [lz_val, ~] = LZc_baseline_multishuffle( ...
            X_seg, seed + ch + w, nShuffles);

        temporal_LZc(ch, w) = lz_val;

    end
end

end