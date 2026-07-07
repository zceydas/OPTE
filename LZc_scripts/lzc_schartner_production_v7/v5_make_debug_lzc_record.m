function rec = v5_make_debug_lzc_record(lzc, window_index, start_sec, end_sec, debug_mode)

rec = struct();
rec.Window = double(window_index);
rec.WindowStartSec = double(start_sec);
rec.WindowEndSec = double(end_sec);
rec.LZc = double(lzc.LZc);
rec.LZcN = double(lzc.LZcN);
rec.RawLZc = double(lzc.RawLZc);
rec.BinaryShuffleMeanRawLZc = double(lzc.BinaryShuffleMeanRawLZc);
rec.BinaryShuffleSDRawLZc = double(lzc.BinaryShuffleSDRawLZc);
rec.BinaryShuffleRawLZcValues = double(lzc.BinaryShuffleRawLZcValues(:));
rec.PhaseRawLZcMean = double(lzc.PhaseRawLZcMean);
rec.PhaseRawLZcSD = double(lzc.PhaseRawLZcSD);
rec.PhaseRawLZcValues = double(lzc.PhaseRawLZcValues(:));
rec.PropOnes = double(lzc.PropOnes);
rec.NTransitions = double(lzc.NTransitions);
rec.StringLength = double(lzc.StringLength);
rec.BinaryStringFirst1000 = string(lzc.BinaryStringFirst1000);

if string(debug_mode) == "full"
    rec.B = lzc.B;
    rec.M = lzc.M;
    rec.TH = lzc.TH;
end

end
