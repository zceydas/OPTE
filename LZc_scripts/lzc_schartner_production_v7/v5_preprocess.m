function Z = v5_preprocess(X)
% Match Schartner python_lzc_py3.py Pre(X):
% Z[i,:] = signal.detrend(X[i,:] - mean(X[i,:]), axis=0)

X = double(X);

if isvector(X)
    X = X(:)';
end

[nChannels, nSamples] = size(X);
Z = zeros(nChannels, nSamples);

for ch = 1:nChannels
    x = X(ch, :);
    Z(ch, :) = detrend(x - mean(x, 'omitnan'));
end

end
