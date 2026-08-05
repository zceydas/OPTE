required_packages <- c("tidyverse", "lme4", "lmerTest", "emmeans", "MuMIn")
missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them in R with:\ninstall.packages(c(\"",
    paste(missing_packages, collapse = "\", \""),
    "\"))\n",
    call. = FALSE
  )
}

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)

# =============================================================================
# OPTE mixed-effects starting point: PSD band changes + channel-wise LZs
#
# Broad goal:
#   1. Test how PSD bands change across baseline, dosing epochs, and eyes.
#   2. Build an analysis-ready PSD + LZs table for later MEQ/subjective models.
#   3. Run exploratory PSD-LZs association models across all bands, not only alpha.
#
# Main PSD model:
#   z_PSD_measure ~ Eyes * Epoch + (1 | Participant)
#
# Baseline-relative PSD model:
#   Delta_from_baseline ~ Eyes * Epoch + (1 | Participant)
#
# Exploratory PSD-LZs model:
#   z_PSD_measure ~ z_LZs + Eyes + Epoch + (1 | Participant)
#
# Notes:
#   - PSD is averaged across channels before modeling.
#   - LZs is averaged across channels and windows within Participant x Eyes x Epoch.
#   - Baseline-relative deltas are computed within Participant x Eyes x Band.
# =============================================================================

# ---- Paths -------------------------------------------------------------------
# This script lives in "<repo>/Mixed effects". The default input paths below are
# relative to the Codex copy of the OPTE outputs. Edit these if you move files.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
if (length(file_arg) > 0) {
  script_file <- tryCatch(
    normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE),
    error = function(e) NA_character_
  )
} else {
  script_file <- tryCatch(
    normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
    error = function(e) NA_character_
  )
}
if (is.na(script_file) || script_file == "") {
  script_dir <- getwd()
} else {
  script_dir <- dirname(script_file)
}
script_dir <- gsub("~+~", " ", script_dir, fixed = TRUE)
script_dir <- normalizePath(script_dir, mustWork = FALSE)
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

lzs_csv <- file.path(
  repo_root,
  "MATLAB/OPTE/data/outputs/resting_preprocessing",
  "LZc_Schartner_ProductionV7_Results/CSV/PUBLIC_long_LZs_by_channel_window.csv"
)

psd_band_csv <- file.path(
  repo_root,
  "MATLAB/OPTE/data/outputs/resting_preprocessing",
  "PSD_Results/ALL_PSD_band_power_by_channel.csv"
)

env_lzs_csv <- Sys.getenv("OPTE_LZS_CSV", unset = "")
env_psd_band_csv <- Sys.getenv("OPTE_PSD_BAND_CSV", unset = "")
if (nzchar(env_lzs_csv)) {
  lzs_csv <- env_lzs_csv
}
if (nzchar(env_psd_band_csv)) {
  psd_band_csv <- env_psd_band_csv
}

output_dir <- file.path(script_dir, "outputs", "PSD_band_LZs_mixed_models")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found:\n", path, call. = FALSE)
  }
}

stop_if_missing(lzs_csv, "LZs CSV")
stop_if_missing(psd_band_csv, "PSD band CSV")

cat("\nUsing PSD band CSV:\n", psd_band_csv, "\n", sep = "")
cat("Using LZs CSV:\n", lzs_csv, "\n", sep = "")
cat("Writing outputs to:\n", output_dir, "\n\n", sep = "")

# ---- Helpers -----------------------------------------------------------------

clean_epoch <- function(x) {
  x <- str_trim(as.character(x))
  x <- if_else(str_detect(str_to_lower(x), "^epoch\\d+$"),
               str_replace(str_to_lower(x), "epoch", "Epoch"),
               x)
  factor(x, levels = paste0("Epoch", 0:4))
}

clean_ids <- function(df) {
  df %>%
    mutate(
      Participant = str_pad(as.character(Participant), width = 3, pad = "0"),
      Session = str_to_lower(str_trim(as.character(Session))),
      Eyes = str_to_upper(str_trim(as.character(Eyes))),
      Epoch = clean_epoch(Epoch)
    )
}

z_within <- function(x) {
  if (sum(is.finite(x)) < 2 || sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

partial_eta2_from_anova <- function(aov_tbl) {
  aov_tbl <- as.data.frame(aov_tbl)
  aov_tbl$Term <- rownames(aov_tbl)
  rownames(aov_tbl) <- NULL

  if (all(c("F value", "NumDF", "DenDF") %in% names(aov_tbl))) {
    aov_tbl$eta2p <- with(aov_tbl, (`F value` * NumDF) /
                            (`F value` * NumDF + DenDF))
  }

  as_tibble(aov_tbl)
}

fit_lmer_safely <- function(data, formula_text) {
  tryCatch(
    lmer(as.formula(formula_text), data = data, REML = TRUE),
    error = function(e) e
  )
}

summarise_model <- function(model, model_name, band, measure, extra = list()) {
  if (inherits(model, "error")) {
    return(tibble(
      Model = model_name,
      Band = band,
      Measure = measure,
      Term = NA_character_,
      Error = model$message
    ))
  }

  aov_tbl <- partial_eta2_from_anova(anova(model))
  r2 <- suppressWarnings(r.squaredGLMM(model))

  aov_tbl %>%
    mutate(
      Model = model_name,
      Band = band,
      Measure = measure,
      R2_marginal = as.numeric(r2[1]),
      R2_conditional = as.numeric(r2[2]),
      Sigma = sigma(model),
      Error = NA_character_,
      !!!extra
    ) %>%
    relocate(Model, Band, Measure, Term)
}

summarise_fixed_effects <- function(model, model_name, band, measure,
                                    extra = list()) {
  if (inherits(model, "error")) {
    return(tibble(
      Model = model_name,
      Band = band,
      Measure = measure,
      Term = NA_character_,
      Error = model$message
    ))
  }

  coef_tbl <- as.data.frame(summary(model)$coefficients)
  coef_tbl$Term <- rownames(coef_tbl)
  rownames(coef_tbl) <- NULL

  as_tibble(coef_tbl) %>%
    mutate(
      Model = model_name,
      Band = band,
      Measure = measure,
      Error = NA_character_,
      !!!extra
    ) %>%
    relocate(Model, Band, Measure, Term)
}

summarise_emmeans <- function(model, model_name, band, measure, extra = list()) {
  if (inherits(model, "error")) {
    return(tibble())
  }

  tryCatch({
    as_tibble(emmeans(model, ~ Eyes * Epoch)) %>%
      mutate(
        Model = model_name,
        Band = band,
        Measure = measure,
        !!!extra
      ) %>%
      relocate(Model, Band, Measure)
  }, error = function(e) tibble())
}

# ---- Read and summarize PSD --------------------------------------------------

psd_raw <- read_csv(psd_band_csv, show_col_types = FALSE,
                    col_types = cols(Participant = col_character()))

required_psd_cols <- c(
  "Participant", "Session", "Eyes", "Epoch", "Channel", "Band",
  "AbsolutePower", "Log10Power", "RelativePower"
)
missing_psd_cols <- setdiff(required_psd_cols, names(psd_raw))
if (length(missing_psd_cols) > 0) {
  stop("PSD CSV is missing columns: ", paste(missing_psd_cols, collapse = ", "))
}

psd_epoch <- psd_raw %>%
  clean_ids() %>%
  mutate(Band = factor(Band, levels = c("delta", "theta", "alpha", "beta", "gamma"))) %>%
  group_by(Participant, Session, Eyes, Epoch, Band) %>%
  summarise(
    AbsolutePower = mean(AbsolutePower, na.rm = TRUE),
    Log10Power = mean(Log10Power, na.rm = TRUE),
    RelativePower = mean(RelativePower, na.rm = TRUE),
    N_Channels_PSD = n_distinct(Channel),
    .groups = "drop"
  ) %>%
  group_by(Participant, Band) %>%
  mutate(
    z_AbsolutePower = z_within(AbsolutePower),
    z_Log10Power = z_within(Log10Power),
    z_RelativePower = z_within(RelativePower)
  ) %>%
  ungroup()

psd_baseline <- psd_epoch %>%
  filter(Session == "baseline", Epoch == "Epoch0") %>%
  select(
    Participant, Eyes, Band,
    Baseline_AbsolutePower = AbsolutePower,
    Baseline_Log10Power = Log10Power,
    Baseline_RelativePower = RelativePower
  )

psd_epoch <- psd_epoch %>%
  left_join(psd_baseline, by = c("Participant", "Eyes", "Band")) %>%
  mutate(
    Delta_AbsolutePower = AbsolutePower - Baseline_AbsolutePower,
    Delta_Log10Power = Log10Power - Baseline_Log10Power,
    Delta_RelativePower = RelativePower - Baseline_RelativePower
  ) %>%
  group_by(Participant, Band) %>%
  mutate(
    z_Delta_AbsolutePower = z_within(Delta_AbsolutePower),
    z_Delta_Log10Power = z_within(Delta_Log10Power),
    z_Delta_RelativePower = z_within(Delta_RelativePower)
  ) %>%
  ungroup()

write_csv(psd_epoch, file.path(output_dir, "analysis_ready_psd_band_by_epoch.csv"))

# ---- Read and summarize LZs --------------------------------------------------

lzs_raw <- read_csv(lzs_csv, show_col_types = FALSE,
                    col_types = cols(Participant = col_character()))

required_lzs_cols <- c(
  "Participant", "Session", "Eyes", "Epoch", "WindowSec", "Channel",
  "LZs", "LZsN"
)
missing_lzs_cols <- setdiff(required_lzs_cols, names(lzs_raw))
if (length(missing_lzs_cols) > 0) {
  stop("LZs CSV is missing columns: ", paste(missing_lzs_cols, collapse = ", "))
}

lzs_epoch <- lzs_raw %>%
  clean_ids() %>%
  group_by(Participant, Session, Eyes, Epoch, WindowSec) %>%
  summarise(
    Mean_LZs = mean(LZs, na.rm = TRUE),
    Mean_LZsN = mean(LZsN, na.rm = TRUE),
    N_Channels_LZs = n_distinct(Channel),
    N_Windows_LZs = n_distinct(Window),
    .groups = "drop"
  ) %>%
  group_by(Participant, WindowSec) %>%
  mutate(
    z_Mean_LZs = z_within(Mean_LZs),
    z_Mean_LZsN = z_within(Mean_LZsN)
  ) %>%
  ungroup()

write_csv(lzs_epoch, file.path(output_dir, "analysis_ready_lzs_by_epoch_window.csv"))

# ---- Merge PSD and LZs -------------------------------------------------------

psd_lzs <- psd_epoch %>%
  left_join(
    lzs_epoch,
    by = c("Participant", "Session", "Eyes", "Epoch")
  ) %>%
  mutate(
    Participant = factor(Participant),
    Eyes = factor(Eyes, levels = c("EC", "EO")),
    Epoch = factor(Epoch, levels = paste0("Epoch", 0:4)),
    Band = factor(Band, levels = c("delta", "theta", "alpha", "beta", "gamma")),
    WindowSec = factor(WindowSec)
  )

write_csv(psd_lzs, file.path(output_dir, "analysis_ready_psd_band_lzs_by_epoch_window.csv"))

cat("\nAnalysis-ready rows:\n")
cat("  PSD epoch rows:     ", nrow(psd_epoch), "\n", sep = "")
cat("  LZs epoch rows:     ", nrow(lzs_epoch), "\n", sep = "")
cat("  PSD + LZs rows:     ", nrow(psd_lzs), "\n", sep = "")
cat("  Participants:       ", n_distinct(psd_epoch$Participant), "\n", sep = "")
cat("  Bands:              ", paste(levels(droplevels(psd_epoch$Band)), collapse = ", "), "\n", sep = "")

# ---- Mixed models: PSD band level -------------------------------------------

level_measures <- c("z_AbsolutePower", "z_Log10Power", "z_RelativePower")
delta_measures <- c("Delta_AbsolutePower", "Delta_Log10Power", "Delta_RelativePower")

level_anova <- list()
level_fixed <- list()
level_emm <- list()

for (band_name in levels(droplevels(psd_epoch$Band))) {
  for (measure in level_measures) {
    model_data <- psd_epoch %>%
      filter(Band == band_name, is.finite(.data[[measure]])) %>%
      mutate(
        Participant = factor(Participant),
        Eyes = factor(Eyes, levels = c("EC", "EO")),
        Epoch = factor(Epoch, levels = paste0("Epoch", 0:4))
      )

    model <- fit_lmer_safely(model_data, paste0(measure, " ~ Eyes * Epoch + (1 | Participant)"))

    level_anova[[paste(band_name, measure)]] <- summarise_model(
      model, "PSD level: Eyes * Epoch", band_name, measure,
      extra = list(N = nrow(model_data), N_Participants = n_distinct(model_data$Participant))
    )
    level_fixed[[paste(band_name, measure)]] <- summarise_fixed_effects(
      model, "PSD level: Eyes * Epoch", band_name, measure,
      extra = list(N = nrow(model_data), N_Participants = n_distinct(model_data$Participant))
    )
    level_emm[[paste(band_name, measure)]] <- summarise_emmeans(
      model, "PSD level: Eyes * Epoch", band_name, measure
    )
  }
}

level_anova <- bind_rows(level_anova)
level_fixed <- bind_rows(level_fixed)
level_emm <- bind_rows(level_emm)

write_csv(level_anova, file.path(output_dir, "model_results_psd_level_anova.csv"))
write_csv(level_fixed, file.path(output_dir, "model_results_psd_level_fixed_effects.csv"))
write_csv(level_emm, file.path(output_dir, "model_results_psd_level_emmeans.csv"))

# ---- Mixed models: PSD baseline-relative change ------------------------------

delta_anova <- list()
delta_fixed <- list()
delta_emm <- list()

for (band_name in levels(droplevels(psd_epoch$Band))) {
  for (measure in delta_measures) {
    model_data <- psd_epoch %>%
      filter(
        Band == band_name,
        Epoch != "Epoch0",
        is.finite(.data[[measure]])
      ) %>%
      mutate(
        Participant = factor(Participant),
        Eyes = factor(Eyes, levels = c("EC", "EO")),
        Epoch = factor(Epoch, levels = paste0("Epoch", 1:4))
      )

    model <- fit_lmer_safely(model_data, paste0(measure, " ~ Eyes * Epoch + (1 | Participant)"))

    delta_anova[[paste(band_name, measure)]] <- summarise_model(
      model, "PSD delta from baseline: Eyes * Epoch", band_name, measure,
      extra = list(N = nrow(model_data), N_Participants = n_distinct(model_data$Participant))
    )
    delta_fixed[[paste(band_name, measure)]] <- summarise_fixed_effects(
      model, "PSD delta from baseline: Eyes * Epoch", band_name, measure,
      extra = list(N = nrow(model_data), N_Participants = n_distinct(model_data$Participant))
    )
    delta_emm[[paste(band_name, measure)]] <- summarise_emmeans(
      model, "PSD delta from baseline: Eyes * Epoch", band_name, measure
    )
  }
}

delta_anova <- bind_rows(delta_anova)
delta_fixed <- bind_rows(delta_fixed)
delta_emm <- bind_rows(delta_emm)

write_csv(delta_anova, file.path(output_dir, "model_results_psd_delta_anova.csv"))
write_csv(delta_fixed, file.path(output_dir, "model_results_psd_delta_fixed_effects.csv"))
write_csv(delta_emm, file.path(output_dir, "model_results_psd_delta_emmeans.csv"))

# ---- Exploratory models: PSD-LZs association ---------------------------------
# This is not yet the subjective-effect model. It just asks whether epochs with
# higher/lower LZs also show higher/lower PSD band power, controlling for eyes
# and epoch. We run it separately for 2s and 10s LZs windows.

association_anova <- list()
association_fixed <- list()

for (band_name in levels(droplevels(psd_lzs$Band))) {
  for (measure in level_measures) {
    for (window_sec in levels(droplevels(psd_lzs$WindowSec))) {
      model_data <- psd_lzs %>%
        filter(
          Band == band_name,
          WindowSec == window_sec,
          is.finite(.data[[measure]]),
          is.finite(z_Mean_LZs)
        ) %>%
        mutate(
          Participant = factor(Participant),
          Eyes = factor(Eyes, levels = c("EC", "EO")),
          Epoch = factor(Epoch, levels = paste0("Epoch", 0:4))
        )

      model <- fit_lmer_safely(
        model_data,
        paste0(measure, " ~ z_Mean_LZs + Eyes + Epoch + (1 | Participant)")
      )

      key <- paste(band_name, measure, window_sec)
      association_anova[[key]] <- summarise_model(
        model, "PSD predicted by LZs + Eyes + Epoch", band_name, measure,
        extra = list(
          WindowSec = window_sec,
          N = nrow(model_data),
          N_Participants = n_distinct(model_data$Participant)
        )
      )
      association_fixed[[key]] <- summarise_fixed_effects(
        model, "PSD predicted by LZs + Eyes + Epoch", band_name, measure,
        extra = list(
          WindowSec = window_sec,
          N = nrow(model_data),
          N_Participants = n_distinct(model_data$Participant)
        )
      )
    }
  }
}

association_anova <- bind_rows(association_anova)
association_fixed <- bind_rows(association_fixed)

write_csv(association_anova, file.path(output_dir, "model_results_psd_lzs_association_anova.csv"))
write_csv(association_fixed, file.path(output_dir, "model_results_psd_lzs_association_fixed_effects.csv"))

# ---- Compact console summary -------------------------------------------------

cat("\nWrote analysis-ready CSVs and model summaries.\n")
cat("Key files:\n")
cat("  analysis_ready_psd_band_by_epoch.csv\n")
cat("  analysis_ready_lzs_by_epoch_window.csv\n")
cat("  analysis_ready_psd_band_lzs_by_epoch_window.csv\n")
cat("  model_results_psd_level_anova.csv\n")
cat("  model_results_psd_delta_anova.csv\n")
cat("  model_results_psd_lzs_association_fixed_effects.csv\n\n")

cat("Suggested first read:\n")
cat("  1. model_results_psd_level_anova.csv: band differences by Eyes/Epoch\n")
cat("  2. model_results_psd_delta_anova.csv: dosing change from baseline\n")
cat("  3. model_results_psd_lzs_association_fixed_effects.csv: PSD-LZs links\n")
