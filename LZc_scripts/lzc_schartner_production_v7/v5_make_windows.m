function windows = v5_make_windows(nSamples, window_samples, overlap_samples, srate)

step_samples = window_samples - overlap_samples;

if step_samples <= 0
    error('overlap_samples must be less than window_samples.');
end

starts = 1:step_samples:(nSamples - window_samples + 1);
nWindows = numel(starts);

Window = double((1:nWindows)');
StartSample = double(starts(:));
EndSample = double(starts(:) + window_samples - 1);
StartSec = double((StartSample - 1) ./ srate);
EndSec = double(EndSample ./ srate);

windows = table(Window, StartSample, EndSample, StartSec, EndSec);

end
