function [lz_norm, shuffle_mean, shuffle_sd, shuffle_values] = v5_binary_shuffle_norm(raw_lz, binary_string, seed, nShuffles)

if nargin < 3
    seed = [];
end
if nargin < 4 || isempty(nShuffles)
    nShuffles = 1;
end

binary_string = char(binary_string(:)');
shuffle_values = nan(nShuffles, 1);

for s = 1:nShuffles

    if ~isempty(seed)
        rng(seed + s);
    end

    shuffled = binary_string(randperm(numel(binary_string)));
    shuffle_values(s) = double(v5_cpr(shuffled));

end

shuffle_mean = double(mean(shuffle_values, 'omitnan'));
shuffle_sd = double(std(shuffle_values, 'omitnan'));

if isfinite(shuffle_mean) && shuffle_mean ~= 0
    lz_norm = double(raw_lz / shuffle_mean);
else
    lz_norm = NaN;
end

end
