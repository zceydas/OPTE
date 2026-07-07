# LZc Paper-Method Pipeline

This package contains:

- `lzc_paper_method_pipeline.m`

It uses your existing helper functions and follows the paper-style procedure more closely:

1. Non-overlapping 2-second segments.
2. Ordinary LZc/LZs from each segment.
3. Phase-randomized surrogate segment data.
4. Recomputed LZc/LZs on surrogate data.
5. `LZcN = LZc / mean(phase-surrogate LZc)`.
6. File-level summaries averaged across segments.

Important settings near the top:

```matlab
segment_sec = 2;
nBinaryShuffles = 10;
nPhaseSurrogates = 10;

do_LZc_random_channel_picks = true;
nChannelPicks = 30;
nChannelsPerPick = 10;

do_LZs_by_channel = true;
```

For a quick single-file debug run, set:

```matlab
target_participant = "005";
target_session = "baseline";
target_eyes = "EC";
target_epoch = "Epoch0";
nPhaseSurrogates = 1;
```

Outputs are saved under:

```text
LZc_PaperMethod_Results/
```

Main output file:

```text
ALL_paper_method_file_summary.csv
```
## Path Configuration

This script uses the same public path pattern as the other upload scripts:

```matlab
USE_HARDCODED_PATHS = false;  % default: choose folders interactively
```

For batch or reproducible runs, set `USE_HARDCODED_PATHS = true` and fill in `HARD_CODED_EEGLAB_PATH`, `HARD_CODED_INPUT_ROOT`, and optionally `HARD_CODED_OUTPUT_ROOT`.

