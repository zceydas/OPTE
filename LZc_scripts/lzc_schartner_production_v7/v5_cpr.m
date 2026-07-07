function count = v5_cpr(binary_string)
% MATLAB equivalent of Schartner python_lzc_py3.py cpr().

binary_string = char(binary_string(:)');

dict = containers.Map('KeyType', 'char', 'ValueType', 'logical');
w = '';

for i = 1:numel(binary_string)

    c = binary_string(i);
    wc = [w c];

    if isKey(dict, wc)
        w = wc;
    else
        dict(wc) = true;
        w = c;
    end

end

count = double(dict.Count);

end
