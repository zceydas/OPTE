function v5_validate_identity(metric_name, reported, numerator, denominator, context)

reported = double(reported);
numerator = double(numerator);
denominator = double(denominator);

if ~isfinite(reported) && (~isfinite(numerator) || ~isfinite(denominator) || denominator == 0)
    return;
end

expected = numerator / denominator;
tol = 1e-10;

if abs(reported - expected) > tol
    error('%s identity failed for %s. Reported %.15g, expected %.15g = %.15g / %.15g.', ...
        metric_name, context, reported, expected, numerator, denominator);
end

end
