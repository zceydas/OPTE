function [lz_val, out] = LZc_baseline_multishuffle(X, seed, nShuffles)

if nargin < 2
    seed = [];
end

if nargin < 3
    nShuffles = 10;
end

X_pre = pre_lzc(X);
[s, M, TH, B] = str_col_lzc(X_pre);

c_orig = double(cpr_lzc(s));

c_shuf_all = nan(nShuffles, 1);

if ~isempty(seed)
    rng(seed);
end

for sh = 1:nShuffles
    perm_idx = randperm(length(s));
    s_shuf = s(perm_idx);
    c_shuf_all(sh) = double(cpr_lzc(s_shuf));
end

c_shuf_mean = mean(c_shuf_all, 'omitnan');
c_shuf_sd = std(c_shuf_all, 'omitnan');

lz_val = double(c_orig) / double(c_shuf_mean);

out = struct();
out.X_pre = X_pre;
out.M = M;
out.TH = TH;
out.B = B;
out.s = s;
out.c_orig = c_orig;
out.c_shuf_all = c_shuf_all;
out.nShuffles = nShuffles;
out.c_shuf_mean = c_shuf_mean;
out.c_shuf_sd = c_shuf_sd;
out.lz_val = lz_val;

end