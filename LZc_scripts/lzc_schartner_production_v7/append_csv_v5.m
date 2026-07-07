function append_csv_v5(T, csv_path)

if isempty(T)
    return;
end

if exist(csv_path, 'file')
    writetable(T, csv_path, 'WriteMode', 'append', 'WriteVariableNames', false);
else
    writetable(T, csv_path);
end

end
