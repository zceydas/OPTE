function x = asnum(x)

if isnumeric(x)
    x = double(x);
elseif iscell(x)
    x = str2double(string(x));
elseif isstring(x) || ischar(x) || iscategorical(x)
    x = str2double(string(x));
else
    try
        x = double(x);
    catch
        x = str2double(string(x));
    end
end

x = x(:);

end
