function [s, M, TH, B] = str_col_lzc(X)

[ro, co] = size(X);

TH = zeros(ro, 1);
M  = zeros(ro, co);
B  = zeros(ro, co);

for i = 1:ro
    M(i, :) = abs(hilbert(X(i, :)));
    TH(i) = mean(M(i, :));
end

chars = repmat('0', 1, ro * co);
idx = 1;

for j = 1:co
    for i = 1:ro
        if M(i, j) > TH(i)
            B(i, j) = 1;
            chars(idx) = '1';
        else
            B(i, j) = 0;
            chars(idx) = '0';
        end
        idx = idx + 1;
    end
end

s = chars;

end