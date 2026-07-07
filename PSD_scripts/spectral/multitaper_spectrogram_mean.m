function [spec_power, spec_times, spec_freqs] = multitaper_spectrogram_mean( ...
    X, srate, window_sec, step_sec, freq_vector, tapers)

[nChannels, nSamples] = size(X);

window_samples = round(window_sec * srate);
step_samples = round(step_sec * srate);

start_samples = 1:step_samples:(nSamples - window_samples + 1);
nWindows = numel(start_samples);

spec_freqs = freq_vector;
nFreqs = numel(spec_freqs);

spec_power_channels = nan(nChannels, nFreqs, nWindows);
spec_times = nan(1, nWindows);

NW = tapers(1);
K = tapers(2);

[dpss_tapers, ~] = dpss(window_samples, NW, K);

for w = 1:nWindows

    idx1 = start_samples(w);
    idx2 = idx1 + window_samples - 1;

    spec_times(w) = ((idx1 + idx2) / 2) / srate;

    for ch = 1:nChannels

        x = double(X(ch, idx1:idx2));
        x = x - mean(x, 'omitnan');

        taper_power = nan(K, nFreqs);

        for k = 1:K

            xt = x(:) .* dpss_tapers(:,k);

            nfft = max(2^nextpow2(window_samples), 2048);
            xdft = fft(xt, nfft);

            f = (0:nfft-1) * (srate / nfft);
            keep = f <= srate/2;

            f = f(keep);
            pxx = abs(xdft(keep)).^2 / (srate * sum(dpss_tapers(:,k).^2));

            taper_power(k,:) = interp1(f, pxx, spec_freqs, 'linear', NaN);

        end

        spec_power_channels(ch,:,w) = mean(taper_power, 1, 'omitnan');

    end
end

spec_power = squeeze(mean(spec_power_channels, 1, 'omitnan'));

end