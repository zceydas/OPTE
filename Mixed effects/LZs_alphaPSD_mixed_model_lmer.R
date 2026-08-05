required_packages <- c("tidyverse", "lme4", "lmerTest", "emmeans", "effectsize", "MuMIn")
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
library(effectsize)
library(MuMIn)

# ══════════════════════════════════════════════════════════════════════════════
#  Research questions:
#    1. Does LZs differ by eye condition and dosing epoch?
#    2. Does alpha PSD differ by eye condition and dosing epoch?
#    3. Are alpha PSD and LZs related after accounting for eyes and epoch?
#
#  This script intentionally mirrors the organization of LZc_mixed_model.R:
#    - Load and prepare data
#    - Fit Eyes × Epoch mixed model
#    - Report R², Type III ANOVA + partial eta squared
#    - Report estimated marginal means
#    - Report eye and epoch contrasts
#    - Fit polynomial trend model
#
#  New additions for the PSD/LZs project:
#    - Run the same style model for LZs and alpha PSD.
#    - Add an alpha PSD ~ LZs + Eyes + Epoch model.
#    - Add a baseline-change model during dosing:
#        Delta alpha PSD ~ Delta LZs + Eyes + Epoch
#
#  Input file:
#    outputs/PSD_band_LZs_mixed_models/analysis_ready_psd_band_lzs_by_epoch_window.csv
# ══════════════════════════════════════════════════════════════════════════════

# ── 0. Paths ─────────────────────────────────────────────────────────────────

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
input_csv <- file.path(output_dir, "analysis_ready_psd_band_lzs_by_epoch_window.csv")

if (!file.exists(input_csv)) {
  stop(
    "Input CSV not found:\n",
    input_csv,
    "\nRun build_PSD_LZs_mixed_dataset_baseR.R first.",
    call. = FALSE
  )
}

# ── 1. Load & prepare ────────────────────────────────────────────────────────

df_all <- read.csv(input_csv, stringsAsFactors = FALSE,
                   colClasses = c(Participant = "character")) %>%
  filter(Band == "alpha") %>%
  mutate(
    Participant = factor(Participant),
    Eyes = factor(Eyes, levels = c("EC", "EO")),
    Epoch = factor(Epoch, levels = paste0("Epoch", 0:4)),
    WindowSec = as.integer(WindowSec)
  )

cat("\n══════════════════════════════════════════════════════\n")
cat(" Analysis-ready alpha PSD + LZs dataset loaded\n")
cat("══════════════════════════════════════════════════════\n")
cat("Input CSV:\n", input_csv, "\n", sep = "")
cat(sprintf("N obs = %d across %d participants\n",
            nrow(df_all), n_distinct(df_all$Participant)))
cat("Windows available: ", paste(sort(unique(df_all$WindowSec)), collapse = ", "), "\n", sep = "")

# ── 2. Helpers copied from the reference-analysis structure ────────────────────────────────

add_cohens_d <- function(contrast_obj, model) {
  s <- sigma(model)
  tbl <- as.data.frame(contrast_obj)
  tbl$d <- tbl$estimate / s
  tbl
}

add_eta2p <- function(aov_table) {
  aov_table$eta2p <- with(aov_table, (`F value` * NumDF) / (`F value` * NumDF + DenDF))
  aov_table
}

print_r2_anova_summary <- function(model, model_label) {
  r2_vals <- r.squaredGLMM(model)

  cat("\n══════════════════════════════════════════════════════\n")
  cat(model_label, "\n")
  cat("══════════════════════════════════════════════════════\n")
  cat(sprintf("\n R²m = %.3f  |  R²c = %.3f\n", r2_vals[1], r2_vals[2]))

  aov_tbl <- add_eta2p(anova(model))
  cat("\n Type III ANOVA + partial eta squared\n")
  print(round(aov_tbl, 4))

  cat("\n Model summary\n")
  print(summary(model))
}

print_eye_epoch_followups <- function(model, outcome_label) {
  emm <- emmeans(model, ~ Eyes * Epoch)

  cat("\n══════════════════════════════════════════════════════\n")
  cat(outcome_label, " — Estimated marginal means\n", sep = "")
  cat("══════════════════════════════════════════════════════\n")
  print(emm)

  emm_eyes <- emmeans(model, ~ Eyes)
  eyes_main <- contrast(emm_eyes, method = "pairwise")

  cat("\n══════════════════════════════════════════════════════\n")
  cat(outcome_label, " — 1. Main effect of Eyes + Cohen's d\n", sep = "")
  cat("══════════════════════════════════════════════════════\n")
  print(add_cohens_d(eyes_main, model))

  eyes_by_epoch_raw <- as.data.frame(
    contrast(emm, method = "pairwise", by = "Epoch")
  )
  eyes_by_epoch_raw$p_bonf <- pmin(eyes_by_epoch_raw$p.value * 5, 1)
  eyes_by_epoch_raw$d <- eyes_by_epoch_raw$estimate / sigma(model)

  cat("\n══════════════════════════════════════════════════════\n")
  cat(outcome_label, " — 2. EC vs EO at each Epoch (Bonferroni x 5) + d\n", sep = "")
  cat("══════════════════════════════════════════════════════\n")
  print(eyes_by_epoch_raw[, c("Epoch", "contrast", "estimate", "SE",
                              "df", "t.ratio", "p.value", "p_bonf", "d")])

  epoch_vs_base <- contrast(emm, method = "trt.vs.ctrl", ref = 1,
                            by = "Eyes", adjust = "bonferroni")

  cat("\n══════════════════════════════════════════════════════\n")
  cat(outcome_label, " — 3. Each Epoch vs Baseline, by Eyes + d\n", sep = "")
  cat("══════════════════════════════════════════════════════\n")
  print(add_cohens_d(epoch_vs_base, model))

  emm_epoch <- emmeans(model, ~ Epoch)
  epoch_vs_base_av <- contrast(emm_epoch, method = "trt.vs.ctrl", ref = 1,
                               adjust = "bonferroni")

  cat("\n══════════════════════════════════════════════════════\n")
  cat(outcome_label, " — 4. Each Epoch vs Baseline, averaged across Eyes + d\n", sep = "")
  cat("══════════════════════════════════════════════════════\n")
  print(add_cohens_d(epoch_vs_base_av, model))

  epoch_pairwise <- contrast(emm_epoch, method = "pairwise", adjust = "bonferroni")

  cat("\n══════════════════════════════════════════════════════\n")
  cat(outcome_label, " — 5. All pairwise Epoch comparisons, averaged across Eyes + d\n", sep = "")
  cat("══════════════════════════════════════════════════════\n")
  print(add_cohens_d(epoch_pairwise, model))

  epoch_pairwise_by_eyes <- contrast(emm, method = "pairwise",
                                     by = "Eyes", adjust = "bonferroni")

  cat("\n══════════════════════════════════════════════════════\n")
  cat(outcome_label, " — 6. All pairwise Epoch comparisons, by Eyes + d\n", sep = "")
  cat("══════════════════════════════════════════════════════\n")
  print(add_cohens_d(epoch_pairwise_by_eyes, model))
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
#  Run mixed models for each LZs window.
# ══════════════════════════════════════════════════════════════════════════════

for (window_sec in sort(unique(df_all$WindowSec))) {
  df <- df_all %>%
    filter(WindowSec == window_sec) %>%
    droplevels()

  cat("\n\n##############################################################################\n")
  cat(" LZs window: ", window_sec, " seconds\n", sep = "")
  cat("##############################################################################\n")
  cat(sprintf("N obs = %d across %d participants\n",
              nrow(df), n_distinct(df$Participant)))

  # ── LZs Eyes x Epoch mixed model ───────────────────────────────────────────
  # Mirrors the original LZc mixed model, but uses V7 LZs.
  lzs_model <- lmer(z_Mean_LZs ~ Eyes * Epoch + (1 | Participant),
                    data = df, REML = TRUE)

  print_r2_anova_summary(
    lzs_model,
    paste0("LZs mixed model: z_Mean_LZs ~ Eyes * Epoch, window ", window_sec, "s")
  )
  print_eye_epoch_followups(lzs_model, paste0("LZs window ", window_sec, "s"))

  # ── LZs polynomial trend model ─────────────────────────────────────────────
  # Tests whether the LZs trajectory has linear/quadratic shape.
  df <- df %>% mutate(Epoch_num = as.numeric(Epoch) - 1)

  lzs_poly_model <- lmer(z_Mean_LZs ~ Eyes * poly(Epoch_num, 2) + (1 | Participant),
                         data = df, REML = TRUE)

  aov_lzs_poly <- add_eta2p(anova(lzs_poly_model))
  r2_lzs_poly <- r.squaredGLMM(lzs_poly_model)

  cat("\n══════════════════════════════════════════════════════\n")
  cat("LZs polynomial trend model: z_Mean_LZs ~ Eyes x poly(Epoch, 2)\n")
  cat("══════════════════════════════════════════════════════\n")
  cat(sprintf("\n R²m = %.3f  |  R²c = %.3f\n\n", r2_lzs_poly[1], r2_lzs_poly[2]))
  print(round(aov_lzs_poly, 4))
  cat("\n Fixed-effect coefficients\n")
  print(round(as.data.frame(summary(lzs_poly_model)$coefficients), 4))

  # ── Alpha PSD Eyes x Epoch mixed model ─────────────────────────────────────
  # Same model family, now using alpha log10 PSD.
  psd_model <- lmer(z_Log10Power ~ Eyes * Epoch + (1 | Participant),
                    data = df, REML = TRUE)

  print_r2_anova_summary(
    psd_model,
    paste0("Alpha PSD mixed model: z_Log10Power ~ Eyes * Epoch, window ", window_sec, "s")
  )
  print_eye_epoch_followups(psd_model, paste0("Alpha PSD window ", window_sec, "s"))

  # ── Alpha PSD polynomial trend model ───────────────────────────────────────
  psd_poly_model <- lmer(z_Log10Power ~ Eyes * poly(Epoch_num, 2) + (1 | Participant),
                         data = df, REML = TRUE)

  aov_psd_poly <- add_eta2p(anova(psd_poly_model))
  r2_psd_poly <- r.squaredGLMM(psd_poly_model)

  cat("\n══════════════════════════════════════════════════════\n")
  cat("Alpha PSD polynomial trend model: z_Log10Power ~ Eyes x poly(Epoch, 2)\n")
  cat("══════════════════════════════════════════════════════\n")
  cat(sprintf("\n R²m = %.3f  |  R²c = %.3f\n\n", r2_psd_poly[1], r2_psd_poly[2]))
  print(round(aov_psd_poly, 4))
  cat("\n Fixed-effect coefficients\n")
  print(round(as.data.frame(summary(psd_poly_model)$coefficients), 4))

  # ── Alpha PSD predicted by LZs + Eyes + Epoch ──────────────────────────────
  # This is the newer analysis needed for the PSD/LZs project.
  # It asks whether LZs predicts alpha power after accounting for eyes and time.
  psd_from_lzs_model <- lmer(z_Log10Power ~ z_Mean_LZs + Eyes + Epoch + (1 | Participant),
                             data = df, REML = TRUE)

  print_r2_anova_summary(
    psd_from_lzs_model,
    paste0("Alpha PSD from LZs model, window ", window_sec, "s: z_Log10Power ~ z_Mean_LZs + Eyes + Epoch")
  )
  print_predictor_effect(
    psd_from_lzs_model,
    "z_Mean_LZs",
    paste0("Alpha PSD from LZs model, window ", window_sec, "s — z_Mean_LZs effect")
  )

  # ── Baseline-change model during dosing only ───────────────────────────────
  # This asks whether change in LZs from baseline predicts change in alpha PSD.
  dosing_df <- df %>%
    filter(Session == "dosing", Epoch != "Epoch0") %>%
    mutate(Epoch = factor(Epoch, levels = paste0("Epoch", 1:4))) %>%
    droplevels()

  psd_delta_from_lzs_delta_model <- lmer(
    Delta_Log10Power ~ Delta_Mean_LZs + Eyes + Epoch + (1 | Participant),
    data = dosing_df, REML = TRUE
  )

  print_r2_anova_summary(
    psd_delta_from_lzs_delta_model,
    paste0("Alpha PSD change from LZs change, window ", window_sec, "s")
  )
  print_predictor_effect(
    psd_delta_from_lzs_delta_model,
    "Delta_Mean_LZs",
    paste0("Alpha PSD change model, window ", window_sec, "s — Delta_Mean_LZs effect")
  )
}
