LZc Measurement Notes

This folder contains MATLAB scripts for the Schartner-style Lempel-Ziv complexity workflow. The upload package includes both the main V7 full-channel workflow and a Schartner-comparison random 10-channel paper-method workflow.

The scripts produce four related measures:

- `LZs`: single-channel LZ complexity normalized by binary shuffling.
- `LZsN`: single-channel LZ complexity normalized by phase-randomized surrogate data.
- `LZc`: multichannel LZ complexity normalized by binary shuffling.
- `LZcN`: multichannel LZ complexity normalized by phase-randomized surrogate data.


Why There Are Multiple Measures

The single-channel measures, `LZs` and `LZsN`, preserve topographic information. They are useful for heatmaps, topoplots, channel-level summaries, and questions about where complexity differs over the scalp.

The multichannel measures, `LZc` and `LZcN`, combine multiple channels into one binary sequence before estimating complexity. In the main V7 branch this is done using all available channels. In the paper-method branch this is done using repeated random picks of 10 channels to match the Schartner-style comparison.

The binary-shuffle-normalized measures, `LZs` and `LZc`, ask how complex the observed binary sequence is relative to shuffled versions of that sequence.

The phase-normalized measures, `LZsN` and `LZcN`, ask how complex the observed signal is relative to phase-randomized surrogate data. These values are ratios and can be greater than 1 when the observed data have higher complexity than the surrogate denominator.


Random 10-Channel Schartner-Comparison Branch

The `paper_method_random10/` folder contains the repeated random 10-channel workflow. Its default settings are intended to mirror the Schartner et al. 2017 channel selection strategy:

```matlab
do_LZc_random_channel_picks = true;
nChannelPicks = 30;
nChannelsPerPick = 10;
segment_sec = 2;
```

Use this branch when methodological comparability with Schartner et al. 2017 is the priority.

## Windowed and Full-Recording Outputs

The main V7 pipeline calculates 2 s and 10 s windowed outputs. These are useful for comparable fixed-duration estimates and for averaging across windows.

The full-recording supplement calculates one set of values across each complete recording. These outputs are useful as a robustness check and for analyses where the whole resting recording is the unit of interest.

Recommended Interpretation

Use `LZs`/`LZsN` for channel-level and scalp-distribution questions.

Use random 10-channel `LZc`/`LZcN` for Schartner-method comparability.

Use all-channel `LZc`/`LZcN` as a full-channel global complexity extension or robustness check.

Report whether values are windowed or full-recording, whether multichannel values used all channels or random 10-channel subsets, and whether normalization used binary shuffling or phase-randomized surrogates.


WARNING!!!!!!!!!!!!

All LZc scripts are going to take a long time to complete. Running all participants, including all epochs (Baseline, dosing, and all follow-ups for every resting EEG recording) and eye conditions, will take anywhere between 24-48 hours. Redundancies are in place in case the program crashes to ensure your data is saved and any ones already completed are skipped when reinitiating the program.