function v5_posthoc_validate_csvs(csv_root)
% Optional: run this after a pipeline run to validate saved CSVs.

if nargin < 1 || isempty(csv_root)
    csv_root = uigetdir(pwd, 'Select LZc_SchartnerV5_SummarySafe_Results/CSV folder');
end

file_summary_path = fullfile(csv_root, 'ALL_file_summary.csv');
window_summary_path = fullfile(csv_root, 'ALL_window_summary.csv');

F = readtable(file_summary_path, 'TextType', 'string');
W = readtable(window_summary_path, 'TextType', 'string');

for i = 1:height(F)

    subset = W( ...
        W.Participant == F.Participant(i) & ...
        W.Session == F.Session(i) & ...
        W.Eyes == F.Eyes(i) & ...
        W.Epoch == F.Epoch(i) & ...
        W.File == F.File(i) & ...
        W.WindowSec == F.WindowSec(i), :);

    if isempty(subset)
        error('No window rows found for file summary row %d.', i);
    end

    expected_lzcn = mean(asnum(subset.RawLZc) ./ asnum(subset.PhaseRawLZcMean), 'omitnan');
    actual_lzcn = double(F.Mean_LZcN(i));

    if abs(actual_lzcn - expected_lzcn) > 1e-10
        error('Mean_LZcN mismatch at row %d. Actual %.15g, expected %.15g.', i, actual_lzcn, expected_lzcn);
    end

end

fprintf('Posthoc CSV validation passed for %d file-summary rows.\n', height(F));

end
