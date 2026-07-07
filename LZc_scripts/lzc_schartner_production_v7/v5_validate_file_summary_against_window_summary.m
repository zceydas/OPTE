function v5_validate_file_summary_against_window_summary(file_summary, window_summary)

assert_close(double(file_summary.Mean_LZc), mean(asnum(window_summary.LZc), 'omitnan'), 'File Mean_LZc mismatch');
assert_close(double(file_summary.Mean_LZcN), mean(asnum(window_summary.LZcN), 'omitnan'), 'File Mean_LZcN mismatch');
assert_close(double(file_summary.Mean_RawLZc), mean(asnum(window_summary.RawLZc), 'omitnan'), 'File Mean_RawLZc mismatch');
assert_close(double(file_summary.Mean_PhaseRawLZc), mean(asnum(window_summary.PhaseRawLZcMean), 'omitnan'), 'File Mean_PhaseRawLZc mismatch');

% This catches exactly the previous bug: file-level LZcN should be the mean
% of per-window LZcN values, not rounded/coerced/inferred.
expected_lzcn = mean(asnum(window_summary.RawLZc) ./ asnum(window_summary.PhaseRawLZcMean), 'omitnan');
assert_close(double(file_summary.Mean_LZcN), expected_lzcn, 'File Mean_LZcN raw-ratio mismatch');

end

function assert_close(a, b, msg)
tol = 1e-10;
if abs(a - b) > tol
    error('%s. %.15g vs %.15g', msg, a, b);
end
end
