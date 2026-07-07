function [psd_abs, freqs] = welch_psd_by_channel( ...
    X, srate, window_sec, overlap_fraction, nfft_min, freq_min, freq_max)

[nChannels, ~] = size(X);

window_samples = round(window_sec * srate);
window_samples = max(window_samples, 8);

noverlap = round(window_samples * overlap_fraction);
nfft = max(nfft_min, 2^nextpow2(window_samples));

psd_abs = [];
freqs = [];

for ch = 1:nChannels

    x = double(X(ch,:));
    x = x - mean(x, 'omitnan');

    [pxx, f] = pwelch(x, hamming(window_samples), noverlap, nfft, srate);

    keep = f >= freq_min & f <= freq_max;

    if ch == 1
        freqs = f(keep)';
        psd_abs = nan(nChannels, numel(freqs));
    end

    psd_abs(ch,:) = pxx(keep)';

end

end