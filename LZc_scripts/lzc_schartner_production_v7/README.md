# LZc Schartner Production V7 Pipeline

Run:

```matlab
run_lzc_schartner_production_v7
```

## Purpose

This production folder contains the main V7 LZ workflow and the companion reporting scripts. The upload package also includes a separate `paper_method_random10/` folder for the repeated random 10-channel procedure used to connect the analysis back to Schartner-style methodology.

The V7 workflow preserves the V6 long-format/debug structure and adds:

- explicit `PipelineVersion`
- `DateProcessed`
- duplicated `WindowLengthSec` in public CSVs
- public long-format outputs for collaborators
- optional compact/full debug `.mat` files
- validation checks from V5/V6

## Core Method

The shared Schartner-style preprocessing and binary conversion is:

```text
mean-center + detrend
Hilbert envelope
mean threshold
column-wise binary string
Schartner cpr raw LZ
binary-shuffle normalization for LZs/LZc
raw phase-surrogate normalization for LZsN/LZcN
```

There are two complementary multichannel branches in this upload package:

- `run_lzc_schartner_production_v7.m` computes the main V7 windowed outputs using all available channels for all-channel `LZc`/`LZcN`, while also keeping channel-wise `LZs`/`LZsN`.
- `../paper_method_random10/lzc_paper_method_pipeline.m` runs the Schartner-comparison branch with repeated random picks of 10 channels, using settings such as `nChannelPicks = 30` and `nChannelsPerPick = 10`.

The random 10-channel branch is expected when the goal is methodological comparability with Schartner et al. The all-channel V7 branch is retained as a full-channel robustness/extension analysis.

## Measurement Families

This pipeline keeps four related Lempel-Ziv complexity measurements because they answer slightly different questions.

`LZs` is a single-channel, spatially local measure. Each EEG channel is converted into its own binary sequence, Lempel-Ziv complexity is calculated for that channel, and the raw value is divided by the mean raw LZ from binary-shuffled versions of that same sequence. Use `LZs` when the question is about channel-wise/topographic complexity.

`LZsN` is the phase-normalized single-channel version. The numerator is still the raw single-channel LZ value, but the denominator is the mean raw LZ from phase-randomized surrogate data. This is included to ask whether the observed channel-wise complexity is high or low relative to a surrogate that preserves spectral structure while disrupting phase relationships. Values can be greater than 1 when the observed signal is more complex than its phase-surrogate denominator.

`LZc` is the multichannel version. In the main V7 branch, all included channels are binarized and combined into one multichannel binary sequence before computing LZ complexity. In the Schartner-comparison branch, repeated random 10-channel subsets are used instead. In both cases, the value is normalized by binary-shuffled multichannel sequences.

`LZcN` is the phase-normalized multichannel version. It uses the observed multichannel raw LZ numerator, but divides by the mean multichannel raw LZ from phase-randomized surrogate data. Like `LZsN`, values above 1 are possible because this is a ratio against a surrogate baseline, not a bounded percentage.

The pipeline also keeps both windowed and full-recording estimates. The 2 s and 10 s windows make the result comparable across equal-duration segments and allow window-level summaries. The full-recording supplement avoids chopping the recording and gives one estimate per full file, which can be useful as a robustness check or when the scientific question concerns the whole resting-state recording.

## Main V7 Outputs

```text
LZc_Schartner_ProductionV7_Results/CSV/
  parser_check.csv
  ALL_file_summary.csv
  ALL_window_summary.csv
  ALL_channel_LZs_long.csv
  ALL_channel_LZs_summary.csv
  ALL_all_channel_LZc_long.csv
  ALL_window_size_comparison.csv

  PUBLIC_long_LZs_by_channel_window.csv
  PUBLIC_long_LZc_by_window.csv
  PUBLIC_file_summary.csv
```

## Random 10-Channel Paper-Method Outputs

```text
LZc_PaperMethod_Results/
  ALL_paper_method_file_summary.csv
  per-file segment, pick, and channel tables depending on save settings
```

The paper-method branch uses non-overlapping 2 s segments and repeated random 10-channel picks by default.

## Recommended Full Run Settings

```matlab
target_participant = "";
target_session = "";
target_eyes = "";
target_epoch = "";

nBinaryShuffles = 1;      % main V7 branch
nPhaseSurrogates = 10;

save_debug_mat = false;
debug_mode = "compact";
overwrite_existing = false;
```

For the paper-method random 10-channel branch, use:

```matlab
segment_sec = 2;
do_LZc_random_channel_picks = true;
nChannelPicks = 30;
nChannelsPerPick = 10;
nPhaseSurrogates = 10;
```

## Recommended Single-File Test Settings

```matlab
target_participant = "005";
target_session = "baseline";
target_eyes = "EC";
target_epoch = "Epoch0";

nPhaseSurrogates = 1;
save_debug_mat = true;
debug_mode = "compact";
```

## Public CSV Meanings

### `PUBLIC_long_LZs_by_channel_window.csv`

Each row = one participant/session/eyes/epoch/window/channel.

- `LZs`: single-channel binary-shuffle-normalized LZ
- `LZsN`: single-channel raw LZ normalized by phase-surrogate raw LZ
- `RawLZs`: raw Schartner cpr count
- `BinaryShuffleMeanRawLZs`: denominator for LZs
- `PhaseRawLZsMean`: denominator for LZsN

### `PUBLIC_long_LZc_by_window.csv`

Each row = one participant/session/eyes/epoch/window.

- `LZc`: all-channel binary-shuffle-normalized LZ
- `LZcN`: all-channel raw LZ normalized by phase-surrogate raw LZ
- `RawLZc`: raw Schartner cpr count
- `BinaryShuffleMeanRawLZc`: denominator for LZc
- `PhaseRawLZcMean`: denominator for LZcN
- `N_Channels`: all channels used
