function [raw_lz, B, M, TH, binary_string] = v5_binary_string(X)
% Create Schartner-style binary string and raw LZ.

X = v5_preprocess(X);

[nChannels, nSamples] = size(X);

M = zeros(nChannels, nSamples);
TH = zeros(nChannels, 1);
B = false(nChannels, nSamples);

for ch = 1:nChannels
    M(ch, :) = abs(hilbert(X(ch, :)));
    TH(ch) = mean(M(ch, :), 'omitnan');
    B(ch, :) = M(ch, :) > TH(ch);
end

binary_string = reshape(B, [], 1)';
binary_string = char(binary_string + '0');

raw_lz = double(v5_cpr(binary_string));

end
