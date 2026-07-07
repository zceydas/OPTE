function Z = pre_lzc(X)

[ro, co] = size(X);
Z = zeros(ro, co);

for i = 1:ro
    xi = X(i, :) - mean(X(i, :));
    Z(i, :) = detrend(xi);
end

end