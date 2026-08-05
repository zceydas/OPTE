required_packages <- c("tidyverse", "lme4", "lmerTest", "emmeans", "MuMIn")
missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
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

# ══════════════════════════════════════════════════════════════════════════════
#  Research question:
#  Do within-session LZs and alpha PSD predict concurrent MEQ4 scores?
#
#  This script intentionally mirrors the organization of the original
#  LZc_predicts_MEQ4.R script:
#    Model 1: simple association, no time covariate
#    Model 2: time-controlled association, adding Epoch
#    Model 3: interaction model, predictor x Epoch
#    Model comparison: likelihood-ratio tests using ML
#    Model 2 reporting: estimated marginal means and predictor effect size
#
#  Timepoint mapping:
#    MEQ T4 ↔ EEG Epoch1
#    MEQ T6 ↔ EEG Epoch2
#    MEQ T8 ↔ EEG Epoch3
#    MEQ T9 ↔ EEG Epoch4
#
#  Variables:
#    z_Mean_LZs     = V7 Schartner-style LZs, averaged across EC/EO, z-scored
#    z_Log10Power   = alpha-band log10 PSD, averaged across EC/EO, z-scored
#    MEQ4           = concurrent subjective effect score
#
#  Input file:
#    outputs/PSD_band_LZs_mixed_models/alpha_psd_lzs_meq_epoch_avg.csv
#
#  Important:
#    This script requires tidyverse, lme4, lmerTest, emmeans, and MuMIn.
#    If one is missing, install it in R with install.packages().
# ══════════════════════════════════════════════════════════════════════════════

# ── 0. Paths and setup ────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
if (length(file_arg) > 0) {
  script_file <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
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

script_dir <- normalizePath(gsub("~+~", " ", script_dir, fixed = TRUE), mustWork = FALSE)
output_dir <- file.path(script_dir, "outputs", "PSD_band_LZs_mixed_models")
input_csv <- file.path(output_dir, "alpha_psd_lzs_meq_epoch_avg.csv")

if (!file.exists(input_csv)) {
  stop(
    "Input CSV not found:\n",
    input_csv,
    "\nRun alpha_PSD_LZs_MEQ_models_nlme.R first.",
    call. = FALSE
  )
}

# Choose which LZs windows to run. Keep both for comparison with the V7 outputs.
windows_to_run <- c(2, 10)

# ── 1. Load analysis-ready MEQ4 + LZs + alpha PSD data ───────────────────────

df_all <- read.csv(input_csv, stringsAsFactors = FALSE,
                   colClasses = c(Participant = "character")) %>%
  mutate(
    Participant = factor(Participant),
    Epoch = factor(Epoch, levels = paste0("Epoch", 1:4)),
    WindowSec = as.integer(WindowSec)
  ) %>%
  filter(
    !is.na(MEQ4),
    !is.na(z_Mean_LZs),
    !is.na(z_Log10Power)
  )

cat("\n══════════════════════════════════════════════════════\n")
cat(" Analysis-ready dataset loaded\n")
cat("══════════════════════════════════════════════════════\n")
cat("Input CSV:\n", input_csv, "\n", sep = "")
cat(sprintf("N obs = %d across %d participants\n",
            nrow(df_all), n_distinct(df_all$Participant)))
cat("Windows available: ", paste(sort(unique(df_all$WindowSec)), collapse = ", "), "\n", sep = "")

# ── 2. Helpers copied from the reference-analysis structure ────────────────────────────────

add_eta2p <- function(aov_table) {
  aov_table$eta2p <- with(aov_table, (`F value` * NumDF) / (`F value` * NumDF + DenDF))
  aov_table
}

print_model_block <- function(model, model_label) {
  cat("\n══════════════════════════════════════════════════════\n")
  cat(model_label, "\n")
  cat("══════════════════════════════════════════════════════\n")
  cat(sprintf("\n R²m = %.3f  |  R²c = %.3f\n",
              r.squaredGLMM(model)[1], r.squaredGLMM(model)[2]))
  print(summary(model))
  aov <- add_eta2p(anova(model))
  cat(" ANOVA:\n")
  print(round(aov, 4))
}

print_predictor_effect <- function(model, term_name, label) {
  cat("\n══════════════════════════════════════════════════════\n")
  cat(label, "\n")
  cat("══════════════════════════════════════════════════════\n")
  cat(sprintf("  b = %.4f, SE = %.4f, t = %.3f, p = %.4f\n",
              fixef(model)[term_name],
              sqrt(vcov(model)[term_name, term_name]),
              summary(model)$coefficients[term_name, "t value"],
              summary(model)$coefficients[term_name, "Pr(>|t|)"]))
  cat(sprintf("  Cohen's d = %.3f  (b / sigma = %.4f / %.4f)\n",
              fixef(model)[term_name] / sigma(model),
              fixef(model)[term_name], sigma(model)))
}

# ══════════════════════════════════════════════════════════════════════════════
#  Run the sequence separately for 2s and 10s LZs windows.
# ══════════════════════════════════════════════════════════════════════════════

for (window_sec in windows_to_run) {
  df <- df_all %>%
    filter(WindowSec == window_sec) %>%
    droplevels()

  cat("\n\n##############################################################################\n")
  cat(" LZs window: ", window_sec, " seconds\n", sep = "")
  cat("##############################################################################\n")
  print(df %>% select(Participant, Epoch, MEQ_timepoint, MEQ4,
                      Mean_LZs, z_Mean_LZs, Log10Power, z_Log10Power),
        n = Inf)
  cat(sprintf("\nN obs = %d across %d participants\n",
              nrow(df), n_distinct(df$Participant)))

  # ════════════════════════════════════════════════════════════════════════════
  #  LZs MODELS
  # ════════════════════════════════════════════════════════════════════════════

  # MODEL 1 — Simple association: MEQ4 ~ z_Mean_LZs
  # Does LZs predict MEQ4 across the dosing session without time?
  m1_lzs <- lmer(MEQ4 ~ z_Mean_LZs + (1 | Participant), data = df, REML = TRUE)
  print_model_block(
    m1_lzs,
    "LZs Model 1: MEQ4 ~ z_Mean_LZs + (1 | Participant)"
  )

  # MODEL 2 — Time-controlled: MEQ4 ~ z_Mean_LZs + Epoch
  # Does LZs predict MEQ4 above and beyond the shared temporal trend?
  m2_lzs <- lmer(MEQ4 ~ z_Mean_LZs + Epoch + (1 | Participant), data = df, REML = TRUE)
  print_model_block(
    m2_lzs,
    "LZs Model 2: MEQ4 ~ z_Mean_LZs + Epoch + (1 | Participant)"
  )

  # MODEL 3 — Interaction: MEQ4 ~ z_Mean_LZs * Epoch
  # Does the LZs -> MEQ4 relationship vary across timepoints?
  m3_lzs <- lmer(MEQ4 ~ z_Mean_LZs * Epoch + (1 | Participant), data = df, REML = TRUE)
  print_model_block(
    m3_lzs,
    "LZs Model 3: MEQ4 ~ z_Mean_LZs * Epoch + (1 | Participant)"
  )

  # Model comparison for LZs models, using ML for LRT.
  cat("\n══════════════════════════════════════════════════════\n")
  cat("LZs model comparison (likelihood ratio tests, ML)\n")
  cat("══════════════════════════════════════════════════════\n")
  m1_lzs_ml <- lmer(MEQ4 ~ z_Mean_LZs + (1 | Participant), data = df, REML = FALSE)
  m2_lzs_ml <- lmer(MEQ4 ~ z_Mean_LZs + Epoch + (1 | Participant), data = df, REML = FALSE)
  m3_lzs_ml <- lmer(MEQ4 ~ z_Mean_LZs * Epoch + (1 | Participant), data = df, REML = FALSE)
  print(anova(m1_lzs_ml, m2_lzs_ml, m3_lzs_ml))

  # Best LZs model report, matching reference analysis script: Model 2.
  cat("\n══════════════════════════════════════════════════════\n")
  cat("LZs Model 2 — estimated marginal means of MEQ4 by Epoch\n")
  cat("══════════════════════════════════════════════════════\n")
  emm_lzs_m2 <- emmeans(m2_lzs, ~ Epoch)
  print(emm_lzs_m2)
  print_predictor_effect(m2_lzs, "z_Mean_LZs", "LZs Model 2 — z_Mean_LZs effect")

  # ════════════════════════════════════════════════════════════════════════════
  #  ALPHA PSD MODELS
  # ════════════════════════════════════════════════════════════════════════════

  # MODEL 1 — Simple association: MEQ4 ~ z_Log10Power
  # Does alpha PSD predict MEQ4 across the dosing session without time?
  m1_psd <- lmer(MEQ4 ~ z_Log10Power + (1 | Participant), data = df, REML = TRUE)
  print_model_block(
    m1_psd,
    "Alpha PSD Model 1: MEQ4 ~ z_Log10Power + (1 | Participant)"
  )

  # MODEL 2 — Time-controlled: MEQ4 ~ z_Log10Power + Epoch
  # Does alpha PSD predict MEQ4 above and beyond the shared temporal trend?
  m2_psd <- lmer(MEQ4 ~ z_Log10Power + Epoch + (1 | Participant), data = df, REML = TRUE)
  print_model_block(
    m2_psd,
    "Alpha PSD Model 2: MEQ4 ~ z_Log10Power + Epoch + (1 | Participant)"
  )

  # MODEL 3 — Interaction: MEQ4 ~ z_Log10Power * Epoch
  # Does the alpha PSD -> MEQ4 relationship vary across timepoints?
  m3_psd <- lmer(MEQ4 ~ z_Log10Power * Epoch + (1 | Participant), data = df, REML = TRUE)
  print_model_block(
    m3_psd,
    "Alpha PSD Model 3: MEQ4 ~ z_Log10Power * Epoch + (1 | Participant)"
  )

  # Model comparison for alpha PSD models, using ML for LRT.
  cat("\n══════════════════════════════════════════════════════\n")
  cat("Alpha PSD model comparison (likelihood ratio tests, ML)\n")
  cat("══════════════════════════════════════════════════════\n")
  m1_psd_ml <- lmer(MEQ4 ~ z_Log10Power + (1 | Participant), data = df, REML = FALSE)
  m2_psd_ml <- lmer(MEQ4 ~ z_Log10Power + Epoch + (1 | Participant), data = df, REML = FALSE)
  m3_psd_ml <- lmer(MEQ4 ~ z_Log10Power * Epoch + (1 | Participant), data = df, REML = FALSE)
  print(anova(m1_psd_ml, m2_psd_ml, m3_psd_ml))

  # Best alpha PSD model report: Model 2.
  cat("\n══════════════════════════════════════════════════════\n")
  cat("Alpha PSD Model 2 — estimated marginal means of MEQ4 by Epoch\n")
  cat("══════════════════════════════════════════════════════\n")
  emm_psd_m2 <- emmeans(m2_psd, ~ Epoch)
  print(emm_psd_m2)
  print_predictor_effect(m2_psd, "z_Log10Power", "Alpha PSD Model 2 — z_Log10Power effect")

  # ════════════════════════════════════════════════════════════════════════════
  #  COMBINED LZs + ALPHA PSD MODELS
  # ════════════════════════════════════════════════════════════════════════════

  # MODEL 1 — Simple combined association
  # Do LZs and alpha PSD explain MEQ4 together without time?
  m1_combined <- lmer(
    MEQ4 ~ z_Mean_LZs + z_Log10Power + (1 | Participant),
    data = df, REML = TRUE
  )
  print_model_block(
    m1_combined,
    "Combined Model 1: MEQ4 ~ z_Mean_LZs + z_Log10Power + (1 | Participant)"
  )

  # MODEL 2 — Time-controlled combined association
  # Do LZs and alpha PSD explain MEQ4 above and beyond Epoch?
  m2_combined <- lmer(
    MEQ4 ~ z_Mean_LZs + z_Log10Power + Epoch + (1 | Participant),
    data = df, REML = TRUE
  )
  print_model_block(
    m2_combined,
    "Combined Model 2: MEQ4 ~ z_Mean_LZs + z_Log10Power + Epoch + (1 | Participant)"
  )

  # MODEL 3 — Combined interaction model
  # Do the LZs/PSD -> MEQ4 relationships vary across timepoints?
  m3_combined <- lmer(
    MEQ4 ~ (z_Mean_LZs + z_Log10Power) * Epoch + (1 | Participant),
    data = df, REML = TRUE
  )
  print_model_block(
    m3_combined,
    "Combined Model 3: MEQ4 ~ (z_Mean_LZs + z_Log10Power) * Epoch + (1 | Participant)"
  )

  # Model comparison for combined models, using ML for LRT.
  cat("\n══════════════════════════════════════════════════════\n")
  cat("Combined model comparison (likelihood ratio tests, ML)\n")
  cat("══════════════════════════════════════════════════════\n")
  m1_combined_ml <- lmer(
    MEQ4 ~ z_Mean_LZs + z_Log10Power + (1 | Participant),
    data = df, REML = FALSE
  )
  m2_combined_ml <- lmer(
    MEQ4 ~ z_Mean_LZs + z_Log10Power + Epoch + (1 | Participant),
    data = df, REML = FALSE
  )
  m3_combined_ml <- lmer(
    MEQ4 ~ (z_Mean_LZs + z_Log10Power) * Epoch + (1 | Participant),
    data = df, REML = FALSE
  )
  print(anova(m1_combined_ml, m2_combined_ml, m3_combined_ml))

  # Best combined model report: Model 2.
  cat("\n══════════════════════════════════════════════════════\n")
  cat("Combined Model 2 — estimated marginal means of MEQ4 by Epoch\n")
  cat("══════════════════════════════════════════════════════\n")
  emm_combined_m2 <- emmeans(m2_combined, ~ Epoch)
  print(emm_combined_m2)
  print_predictor_effect(m2_combined, "z_Mean_LZs", "Combined Model 2 — z_Mean_LZs effect")
  print_predictor_effect(m2_combined, "z_Log10Power", "Combined Model 2 — z_Log10Power effect")

  # ════════════════════════════════════════════════════════════════════════════
  #  INCREMENTAL COMPARISONS
  # ════════════════════════════════════════════════════════════════════════════

  cat("\n══════════════════════════════════════════════════════\n")
  cat("Incremental test: Does alpha PSD add above LZs + Epoch?\n")
  cat("══════════════════════════════════════════════════════\n")
  print(anova(m2_lzs_ml, m2_combined_ml))

  cat("\n══════════════════════════════════════════════════════\n")
  cat("Incremental test: Does LZs add above alpha PSD + Epoch?\n")
  cat("══════════════════════════════════════════════════════\n")
  print(anova(m2_psd_ml, m2_combined_ml))
}
