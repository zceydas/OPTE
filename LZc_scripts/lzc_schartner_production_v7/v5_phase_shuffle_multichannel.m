function X_surr = v5_phase_shuffle_multichannel(X, seed)

X = double(X);

if isvector(X)
    X = X(:)';
end

[nChannels, nSamples] = size(X);
X_surr = zeros(nChannels, nSamples);

for ch = 1:nChannels
    if isempty(seed)
        ch_seed = [];
    else
        ch_seed = seed + ch;
    end

    X_surr(ch, :) = v5_phase_shuffle_signal(X(ch, :), ch_seed);
end

end
