function nDict = cpr_lzc(s)

d = containers.Map('KeyType', 'char', 'ValueType', 'char');
w = '';

for k = 1:length(s)
    c = s(k);
    wc = [w c];

    if isKey(d, wc)
        w = wc;
    else
        d(wc) = wc;
        w = c;
    end
end

nDict = double(d.Count);

end