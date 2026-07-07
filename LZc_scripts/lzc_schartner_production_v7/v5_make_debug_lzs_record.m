function rec = v5_make_debug_lzs_record(lzs, window_index, start_sec, end_sec, channel_index, channel_label, debug_mode)

rec = struct();
rec.Window = double(window_index);
rec.WindowStartSec = double(start_sec);
rec.WindowEndSec = double(end_sec);
rec.ChannelIndex = double(channel_index);
rec.Channel = string(channel_label);
rec.LZs = double(lzs.LZs);
rec.LZsN = double(lzs.LZsN);
rec.RawLZs = double(lzs.RawLZs);
rec.BinaryShuffleMeanRawLZs = double(lzs.BinaryShuffleMeanRawLZs);
rec.BinaryShuffleSDRawLZs = double(lzs.BinaryShuffleSDRawLZs);
rec.BinaryShuffleRawLZsValues = double(lzs.BinaryShuffleRawLZsValues(:));
rec.PhaseRawLZsMean = double(lzs.PhaseRawLZsMean);
rec.PhaseRawLZsSD = double(lzs.PhaseRawLZsSD);
rec.PhaseRawLZsValues = double(lzs.PhaseRawLZsValues(:));
rec.Threshold = double(lzs.Threshold);
rec.PropOnes = double(lzs.PropOnes);
rec.NTransitions = double(lzs.NTransitions);
rec.StringLength = double(lzs.StringLength);
rec.BinaryStringFirst1000 = string(lzs.BinaryStringFirst1000);

if string(debug_mode) == "full"
    rec.B = lzs.B;
    rec.M = lzs.M;
    rec.TH = lzs.TH;
end

end
