# PSD Power Measures in the MATLAB Scripts

The main script `spectral_analysis.m` computes three power scales from the same Welch PSD output.

## Script Formulas

Absolute power is the PSD returned by `pwelch` in `welch_psd_by_channel.m`.

```matlab
[pxx, f] = pwelch(x, hamming(window_samples), noverlap, nfft, srate);
AbsolutePower = pxx;
```

Log10 power is calculated in `spectral_analysis.m` as:

```matlab
psd_log = log10(psd_abs + eps);
```

Relative power is calculated in `spectral_analysis.m` as:

```matlab
total_power = trapz(freqs, psd_abs, 2);
psd_rel = psd_abs ./ total_power;
```

The total power denominator is computed separately for each channel across the retained frequency range.

## Interpretation

`AbsolutePower` keeps the original PSD magnitude and is useful when raw power differences matter.

`Log10Power` is a transformed version of absolute power. It compresses large values and is often easier to visualize and model statistically.

`RelativePower` normalizes each channel by its total analyzed power, making it useful for comparing the distribution of spectral power across frequencies or bands.

## Where Values Are Written

`psd_to_long_table.m` writes frequency-level rows with `AbsolutePower`, `Log10Power`, and `RelativePower`.

`extract_band_power_table.m` writes band-level rows by averaging each measure within each canonical band.

Group statistics and figures are generated separately for all three measures.
