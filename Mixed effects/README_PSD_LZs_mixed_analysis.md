# OPTE PSD + LZs Mixed-Analysis Starting Point

This folder contains R scripts for moving from cleaned long-format EEG outputs into broad PSD band summaries and mixed-effects models.

The current broad question is:

> How do delta, theta, alpha, beta, and gamma power change between baseline, dosing epochs, and eye conditions?

The alpha-suppression question can be tested later as a narrower follow-up.

## Scripts

1. `build_PSD_LZs_mixed_dataset_baseR.R`

   Builds analysis-ready CSVs from:

   - `PUBLIC_long_LZs_by_channel_window.csv`
   - `ALL_PSD_band_power_by_channel.csv`

   This script only uses base R. Run it first.

2. `summarize_PSD_band_changes_baseR.R`

   Builds descriptive tables and figures for all PSD bands across:

   - baseline Epoch0
   - dosing Epoch1-Epoch4
   - eyes closed and eyes open

   This also only uses base R, so it can run before installing mixed-model packages.

3. `PSD_band_LZs_mixed_models.R`

   Runs the actual mixed-effects models once these packages are installed:

   ```r
   install.packages(c("tidyverse", "lme4", "lmerTest", "emmeans", "MuMIn"))
   ```

4. `PSD_band_LZs_mixed_models_nlme.R`

   Runs a mixed-effects version using `nlme`, which is often already installed
   with R. This is the best immediate option when `lme4`, `lmerTest`, and
   `emmeans` are not installed yet.

## Recommended Run Order

From the repository root:

```bash
Rscript "Mixed effects/build_PSD_LZs_mixed_dataset_baseR.R"
Rscript "Mixed effects/summarize_PSD_band_changes_baseR.R"
Rscript "Mixed effects/PSD_band_LZs_mixed_models_nlme.R"
```

Then, after installing the fuller R modeling packages:

```bash
Rscript "Mixed effects/PSD_band_LZs_mixed_models.R"
```

## Main Outputs

The scripts write to:

`Mixed effects/outputs/PSD_band_LZs_mixed_models`

Useful analysis-ready data:

- `analysis_ready_psd_band_by_epoch.csv`
- `analysis_ready_lzs_by_epoch_window.csv`
- `analysis_ready_psd_band_lzs_by_epoch_window.csv`

Useful descriptive summaries:

- `descriptive_psd_band_summary_by_epoch.csv`
- `descriptive_psd_band_delta_summary_by_epoch.csv`

Useful figures:

- `figures/PSD_band_log10_power_by_epoch_and_eyes.png`
- `figures/PSD_band_relative_power_by_epoch_and_eyes.png`
- `figures/PSD_band_delta_log10_power_by_epoch_and_eyes.png`
- `figures/PSD_band_delta_relative_power_by_epoch_and_eyes.png`

Mixed-model outputs available immediately with `nlme`:

- `nlme_model_results_psd_level_anova.csv`
- `nlme_model_results_psd_level_fixed_effects.csv`
- `nlme_model_results_psd_delta_anova.csv`
- `nlme_model_results_psd_delta_fixed_effects.csv`
- `nlme_model_results_psd_lzs_association_anova.csv`
- `nlme_model_results_psd_lzs_association_fixed_effects.csv`

## What The Measures Mean

Absolute power is the raw amount of power in a band. It is intuitive, but it can be dominated by high-power participants or channels.

Log10 power is the log-transformed version of absolute power. This often behaves better statistically because EEG power is usually skewed.

Relative power is the proportion of total power represented by a band. This is useful for asking whether a band becomes more or less dominant within the spectrum, but it is compositional: if one band goes up, another may appear to go down even if its absolute power did not change much.

Delta-from-baseline columns compare each participant's dosing epochs to their own baseline for the same eye condition and band. These are often the clearest first-pass values for dosing effects.

## Mixed-Effects Model Logic

The broad PSD model is:

```text
PSD measure ~ Eyes * Epoch + (1 | Participant)
```

In words, this asks whether PSD differs by eye condition, epoch, or the combination of eye condition and epoch, while accounting for repeated measurements from the same participant.

The baseline-change model is:

```text
PSD change from baseline ~ Eyes * Epoch + (1 | Participant)
```

This asks whether dosing-related change differs across epochs and eye conditions.

The exploratory PSD-LZs model is:

```text
PSD measure ~ LZs + Eyes + Epoch + (1 | Participant)
```

This asks whether participants/epochs with higher LZs also tend to have different PSD values, after accounting for eyes and epoch.

## Alpha PSD, LZs, And MEQ4 Follow-Up

After the broad band models, the next focused script is:

```bash
Rscript "Mixed effects/alpha_PSD_LZs_MEQ_models_nlme.R"
```

This script focuses on alpha and writes:

- `alpha_psd_lzs_by_eye_window.csv`
- `alpha_psd_lzs_epoch_avg_for_meq.csv`
- `alpha_psd_lzs_meq_epoch_avg.csv`
- `alpha_psd_lzs_correlations.csv`
- `alpha_psd_lzs_nlme_fixed_effects.csv`
- `alpha_meq_psd_lzs_nlme_fixed_effects.csv`

The PSD-LZs models preserve eye condition. The MEQ4 models average across eye
condition because MEQ4 is measured once per dosing timepoint rather than once per
eye condition.

The MEQ4 time mapping follows the existing project script:

- T4 = Epoch1
- T6 = Epoch2
- T8 = Epoch3
- T9 = Epoch4

## MEQ4 Time Models With Alpha PSD

To mirror the structure of the original `LZc_predicts_MEQ4.R` script while
adding alpha PSD:

```bash
Rscript "Mixed effects/LZs_alphaPSD_predicts_MEQ4_time_model_nlme.R"
```

This writes:

- `alpha_psd_lzs_meq_model_dataset.csv`
- `alpha_psd_lzs_meq_fixed_effects.csv`
- `alpha_psd_lzs_meq_anova.csv`
- `alpha_psd_lzs_meq_lrt_model_comparisons.csv`
- `alpha_psd_lzs_meq_incremental_lrt.csv`

The model sequence follows the reference analysis script:

- Model 1: simple association, such as `MEQ4 ~ z_Mean_LZs`
- Model 2: time-controlled association, such as `MEQ4 ~ z_Mean_LZs + Epoch`
- Model 3: time-varying association, such as `MEQ4 ~ z_Mean_LZs * Epoch`

The script runs that sequence for LZs alone, alpha PSD alone, and LZs plus
alpha PSD together. It also tests whether alpha PSD adds to `LZs + Epoch`, and
whether LZs adds to `alpha PSD + Epoch`, using likelihood-ratio comparisons.

Interpretation note: `nlme::anova()` is sequential, so predictor order can
matter when LZs, alpha PSD, and time share variance. For poster/report wording,
use the fixed-effect coefficient table and likelihood-ratio comparisons as the
primary outputs.

To make figures from these time models:

```bash
Rscript "Mixed effects/plot_alphaPSD_LZs_MEQ4_time_models_baseR.R"
```

This writes:

- `figures/meq4_simple_scatter_window2s.png`
- `figures/meq4_simple_scatter_window10s.png`
- `figures/meq4_predictor_attenuation_window2s.png`
- `figures/meq4_predictor_attenuation_window10s.png`
- `figures/meq4_epoch_trajectories_window2s.png`
- `figures/meq4_epoch_trajectories_window10s.png`
- `figures/meq4_incremental_lrt_window2s.png`
- `figures/meq4_incremental_lrt_window10s.png`

The predictor attenuation figure is the most direct visual match to the reference analysis
script's logic: it shows simple EEG-MEQ4 associations, then shows how those
effects weaken after adding `Epoch` as a time covariate. The incremental LRT
figure asks whether alpha PSD adds to `LZs + Epoch`, and whether LZs adds to
`alpha PSD + Epoch`.

For poster figures, use the cleaner plotting script:

```bash
Rscript "Mixed effects/make_poster_figures_baseR.R"
```

This writes white-background, black-text figures with larger significance
labels and legends outside the plotting area:

- `figures/poster/poster_lzs_eye_condition_window2s.png`
- `figures/poster/poster_lzs_eye_condition_window10s.png`
- `figures/poster/poster_alpha_psd_eye_condition.png`
- `figures/poster/poster_meq_time_model_trajectory_window2s.png`
- `figures/poster/poster_meq_time_model_trajectory_window10s.png`
- `figures/poster/poster_meq_predictor_attenuation_window2s.png`
- `figures/poster/poster_meq_predictor_attenuation_window10s.png`
- `figures/poster/poster_meq_incremental_lrt_window2s.png`
- `figures/poster/poster_meq_incremental_lrt_window10s.png`

For poster use, prefer these `poster` files over the earlier draft
figures because the labels and legends were adjusted specifically for slide and
poster layout. The MEQ4 trajectory figures include Baseline/Epoch0 for the EEG
variables, while MEQ4 begins at Epoch1/T4 because MEQ4 was not measured at
baseline. These poster figures label the complexity measure as `LZc`; internal
file names and column names may still use `LZs` because that was the analysis
pipeline's bookkeeping name for this LZc computation. Poster x-axis labels use
study time (`Baseline`, `90 min`, `180 min`, `270 min`, `360 min`) instead of
the internal epoch/timepoint names.

There is also a more literal `lmer` version:

```bash
Rscript "Mixed effects/LZs_alphaPSD_predicts_MEQ4_lmer.R"
```

This script is intentionally organized like the original `LZc_predicts_MEQ4.R`:
it explicitly defines `m1`, `m2`, and `m3` model objects for LZs, alpha PSD, and
combined LZs + alpha PSD models. It requires the same modeling ecosystem as the
reference analysis script:

```r
install.packages(c("tidyverse", "lme4", "lmerTest", "emmeans", "MuMIn"))
```

Use this version when you want the code to be easiest to compare line-by-line
with the reference analysis script. Use the `nlme` version when you want the version that runs
with fewer external package requirements.

For the mixed-model side, there is an equivalent teaching script:

```bash
Rscript "Mixed effects/LZs_alphaPSD_mixed_model_lmer.R"
```

This mirrors `LZc_mixed_model.R` and explicitly reports:

- LZs `Eyes * Epoch` mixed models
- alpha PSD `Eyes * Epoch` mixed models
- estimated marginal means with `emmeans`
- eyes and epoch follow-up contrasts
- polynomial trend models
- alpha PSD predicted by LZs while accounting for eyes and epoch
- alpha PSD change predicted by LZs change during dosing

It requires:

```r
install.packages(c("tidyverse", "lme4", "lmerTest", "emmeans", "effectsize", "MuMIn"))
```

This script is intentionally longer and less automated than
`PSD_band_LZs_mixed_models_nlme.R`, because it is meant to be readable beside
the original reference analysis script while preserving the newer PSD/LZs outputs.

## Line Graph Results

For presentation-style line graphs with error bars:

```bash
Rscript "Mixed effects/plot_alpha_psd_lzs_meq_epoch_lines_baseR.R"
```

This writes:

- `alpha_psd_lzs_meq_epoch_line_summary.csv`
- `alpha_psd_lzs_meq_epoch_correlations.csv`
- `alpha_psd_lzs_with_baseline_epoch_line_summary.csv`
- `alpha_psd_lzs_with_baseline_epoch_correlations.csv`
- `figures/alpha_psd_lzs_with_baseline_epoch_lines_window2s.png`
- `figures/alpha_psd_lzs_with_baseline_epoch_lines_window10s.png`
- `figures/alpha_psd_lzs_with_baseline_epoch_correlations_window2s.png`
- `figures/alpha_psd_lzs_with_baseline_epoch_correlations_window10s.png`
- `figures/alpha_psd_lzs_meq_epoch_lines_window2s.png`
- `figures/alpha_psd_lzs_meq_epoch_lines_window10s.png`
- `figures/alpha_psd_lzs_meq_delta_epoch_lines_window2s.png`
- `figures/alpha_psd_lzs_meq_delta_epoch_lines_window10s.png`
- `figures/alpha_psd_lzs_meq_epoch_correlations_window2s.png`
- `figures/alpha_psd_lzs_meq_epoch_correlations_window10s.png`

The trajectory plots show mean +/- 95% CI across dosing epochs. The correlation
plots show the Pearson correlation at each dosing epoch across participants.

All line graphs now include the Baseline Epoch0 position on the x-axis. Alpha
PSD and LZs have real baseline values. MEQ4 is not measured at baseline, so MEQ4
starts at Epoch1/T4 while the baseline position remains visible for context.

For EEG-only eye-condition plots:

```bash
Rscript "Mixed effects/plot_alpha_eye_condition_contrasts_baseR.R"
```

This writes:

- `alpha_eye_condition_line_summary.csv`
- `alpha_eye_condition_paired_stats.csv`
- `alpha_eo_minus_ec_line_summary.csv`
- `alpha_eo_minus_ec_stats.csv`
- `figures/alpha_psd_eye_condition_lines_window2s.png`
- `figures/alpha_psd_eye_condition_lines_window10s.png`
- `figures/lzs_eye_condition_lines_window2s.png`
- `figures/lzs_eye_condition_lines_window10s.png`
- `figures/alpha_psd_eo_minus_ec_lines_window2s.png`
- `figures/alpha_psd_eo_minus_ec_lines_window10s.png`
- `figures/lzs_eo_minus_ec_lines_window2s.png`
- `figures/lzs_eo_minus_ec_lines_window10s.png`

The EC/EO plots use paired tests at each epoch. The EO-minus-EC plots test
whether the contrast differs from zero at each epoch.
