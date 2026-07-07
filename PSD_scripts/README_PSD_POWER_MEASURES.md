PSD Power Measure Notes

This folder contains the MATLAB scripts for the resting EEG power spectral density workflow. The main runner is:

```matlab
spectral/spectral_analysis.m
```

The PSD pipeline saves three related power measures so downstream users can choose the scale that best matches their analysis question:

- `AbsolutePower`
- `Log10Power`
- `RelativePower`


Absolute Power

`AbsolutePower` is the raw Welch power spectral density estimate from `pwelch`, after selecting the configured frequency range. In this pipeline, the default analysis range is 1-45 Hz.

Absolute power is useful when the question is about the amount of power present at a frequency or within a band. It preserves amplitude information, but it can be strongly affected by overall signal scale, impedance/noise differences, skull/scalp factors, and between-participant differences in total EEG power.

Use absolute power when total signal magnitude is meaningful and should not be normalized away.


Log10 Power

`Log10Power` is calculated as:

```matlab
Log10Power = log10(AbsolutePower + eps);
```

The log transform compresses the large dynamic range of PSD values and often makes power distributions more suitable for visualization and statistics. Differences in log10 power are multiplicative on the original scale: a change of 1 log10 unit corresponds to a 10-fold change in absolute power.

Use log10 power when absolute power is highly skewed or when proportional/multiplicative changes are easier to interpret than raw power differences.


Relative Power

`RelativePower` is calculated per channel by dividing absolute power at each frequency by that channel's total power across the analyzed frequency range:

```matlab
total_power = trapz(freqs, AbsolutePower, 2);
RelativePower = AbsolutePower ./ total_power;
```

Relative power describes how much of a channel's total analyzed power is allocated to each frequency or band. It reduces the influence of overall signal amplitude, but it also creates compositional dependencies: if one frequency band takes a larger share, at least one other part of the spectrum must take a smaller share.

Use relative power when the question is about spectral distribution or proportion of total power rather than raw signal magnitude.


Frequency and Band Outputs

`ALL_PSD_by_channel_frequency.csv` stores one row per participant/session/eyes/epoch/file/channel/frequency. It includes `AbsolutePower`, `Log10Power`, and `RelativePower`.

`ALL_PSD_band_power_by_channel.csv` stores one row per participant/session/eyes/epoch/file/channel/frequency band. Band values are averages across the frequencies inside each band for the same three measures.

Default bands are:

- delta: 1-4 Hz
- theta: 4-8 Hz
- alpha: 8-13 Hz
- beta: 13-30 Hz
- gamma: 30-45 Hz


Practical Recommendation

Report which power scale was used. Absolute, log10, and relative power are not interchangeable:

- absolute power asks "how much power is present?"
- log10 power asks "how does power differ on a compressed multiplicative scale?"
- relative power asks "what fraction of the analyzed spectrum is allocated here?"
