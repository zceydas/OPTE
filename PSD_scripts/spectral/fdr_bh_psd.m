function p_fdr = fdr_bh_psd(p)
p = p(:);
p_fdr = nan(size(p));
valid = isfinite(p);
pv = p(valid);
[ps,idx] = sort(pv);
m = numel(ps);
adj = ps .* m ./ (1:m)';
adj = flipud(cummin(flipud(adj)));
adj(adj>1)=1;
tmp = nan(size(pv));
tmp(idx)=adj;
p_fdr(valid)=tmp;
end
