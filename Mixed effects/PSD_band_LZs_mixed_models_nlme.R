# =============================================================================
# OPTE broad PSD + LZs mixed-effects models using nlme.
#
# This script uses nlme, which is commonly installed with R, as a lightweight
# alternative to the lme4/lmerTest/emmeans script.
#
# Broad PSD models:
#   z_PSD_measure ~ Eyes * Epoch + random intercept for Participant
#
# Baseline-relative dosing models:
#   Delta_PSD_measure ~ Eyes * Epoch + random intercept for Participant
#
# Exploratory PSD-LZs association models:
#   z_PSD_measure ~ z_Mean_LZs + Eyes + Epoch + random intercept for Participant
# =============================================================================

required_packages <- c("nlme")
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

library(nlme)

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
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

psd_csv <- file.path(output_dir, "analysis_ready_psd_band_by_epoch.csv")
combined_csv <- file.path(output_dir, "analysis_ready_psd_band_lzs_by_epoch_window.csv")

env_psd_csv <- Sys.getenv("OPTE_PSD_READY_CSV", unset = "")
env_combined_csv <- Sys.getenv("OPTE_PSD_LZS_READY_CSV", unset = "")
if (nzchar(env_psd_csv)) {
  psd_csv <- env_psd_csv
}
if (nzchar(env_combined_csv)) {
  combined_csv <- env_combined_csv
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found:\n", path, "\nRun build_PSD_LZs_mixed_dataset_baseR.R first.", call. = FALSE)
  }
}

stop_if_missing(psd_csv, "Analysis-ready PSD CSV")
stop_if_missing(combined_csv, "Analysis-ready PSD + LZs CSV")

cat("\nUsing PSD CSV:\n", psd_csv, "\n", sep = "")
cat("Using PSD + LZs CSV:\n", combined_csv, "\n", sep = "")
cat("Writing nlme model outputs to:\n", output_dir, "\n\n", sep = "")

psd <- read.csv(psd_csv, stringsAsFactors = FALSE, colClasses = c(Participant = "character"))
combined <- read.csv(combined_csv, stringsAsFactors = FALSE, colClasses = c(Participant = "character"))

prep_factors <- function(df) {
  df$Participant <- factor(sprintf("%03d", as.integer(as.character(df$Participant))))
  df$Eyes <- factor(toupper(as.character(df$Eyes)), levels = c("EC", "EO"))
  df$Epoch <- factor(as.character(df$Epoch), levels = paste0("Epoch", 0:4))
  df$Band <- factor(tolower(as.character(df$Band)), levels = c("delta", "theta", "alpha", "beta", "gamma"))
  if ("WindowSec" %in% names(df)) {
    df$WindowSec <- factor(df$WindowSec)
  }
  df
}

psd <- prep_factors(psd)
combined <- prep_factors(combined)

fit_lme_safely <- function(data, formula_text) {
  tryCatch(
    lme(
      fixed = as.formula(formula_text),
      random = ~1 | Participant,
      data = data,
      method = "REML",
      na.action = na.omit,
      control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)
    ),
    error = function(e) e
  )
}

anova_table <- function(model, model_name, band, measure, extra = list()) {
  if (inherits(model, "error")) {
    out <- data.frame(
      Model = model_name,
      Band = band,
      Measure = measure,
      Term = NA_character_,
      Error = model$message,
      stringsAsFactors = FALSE
    )
    for (nm in names(extra)) out[[nm]] <- extra[[nm]]
    return(out)
  }

  aov <- as.data.frame(anova(model))
  aov$Term <- rownames(aov)
  rownames(aov) <- NULL
  out <- data.frame(
    Model = model_name,
    Band = band,
    Measure = measure,
    Term = aov$Term,
    numDF = aov$numDF,
    denDF = aov$denDF,
    F_value = aov$`F-value`,
    p_value = aov$`p-value`,
    Error = NA_character_,
    stringsAsFactors = FALSE
  )
  out$partial_eta2 <- with(out, (F_value * numDF) / (F_value * numDF + denDF))
  for (nm in names(extra)) out[[nm]] <- extra[[nm]]
  out
}

fixed_table <- function(model, model_name, band, measure, extra = list()) {
  if (inherits(model, "error")) {
    out <- data.frame(
      Model = model_name,
      Band = band,
      Measure = measure,
      Term = NA_character_,
      Error = model$message,
      stringsAsFactors = FALSE
    )
    for (nm in names(extra)) out[[nm]] <- extra[[nm]]
    return(out)
  }

  coefs <- as.data.frame(summary(model)$tTable)
  coefs$Term <- rownames(coefs)
  rownames(coefs) <- NULL
  out <- data.frame(
    Model = model_name,
    Band = band,
    Measure = measure,
    Term = coefs$Term,
    Estimate = coefs$Value,
    Std_Error = coefs$Std.Error,
    DF = coefs$DF,
    t_value = coefs$`t-value`,
    p_value = coefs$`p-value`,
    Error = NA_character_,
    stringsAsFactors = FALSE
  )
  for (nm in names(extra)) out[[nm]] <- extra[[nm]]
  out
}

append_result <- function(results, item) {
  results[[length(results) + 1]] <- item
  results
}

bind_results <- function(results) {
  if (length(results) == 0) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(results, names)))
  padded <- lapply(results, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[all_names]
  })
  do.call(rbind, padded)
}

level_measures <- c("z_AbsolutePower", "z_Log10Power", "z_RelativePower")
delta_measures <- c("Delta_AbsolutePower", "Delta_Log10Power", "Delta_RelativePower")

level_anova <- list()
level_fixed <- list()

for (band_name in levels(droplevels(psd$Band))) {
  for (measure in level_measures) {
    model_data <- psd[psd$Band == band_name & is.finite(psd[[measure]]), ]
    model_data <- droplevels(model_data)
    model <- fit_lme_safely(model_data, paste0(measure, " ~ Eyes * Epoch"))
    extra <- list(N = nrow(model_data), N_Participants = length(unique(model_data$Participant)))
    level_anova <- append_result(level_anova, anova_table(model, "PSD level: Eyes * Epoch", band_name, measure, extra))
    level_fixed <- append_result(level_fixed, fixed_table(model, "PSD level: Eyes * Epoch", band_name, measure, extra))
  }
}

delta_anova <- list()
delta_fixed <- list()

for (band_name in levels(droplevels(psd$Band))) {
  for (measure in delta_measures) {
    model_data <- psd[psd$Band == band_name & psd$Epoch != "Epoch0" & is.finite(psd[[measure]]), ]
    model_data$Epoch <- factor(as.character(model_data$Epoch), levels = paste0("Epoch", 1:4))
    model_data <- droplevels(model_data)
    model <- fit_lme_safely(model_data, paste0(measure, " ~ Eyes * Epoch"))
    extra <- list(N = nrow(model_data), N_Participants = length(unique(model_data$Participant)))
    delta_anova <- append_result(delta_anova, anova_table(model, "PSD delta from baseline: Eyes * Epoch", band_name, measure, extra))
    delta_fixed <- append_result(delta_fixed, fixed_table(model, "PSD delta from baseline: Eyes * Epoch", band_name, measure, extra))
  }
}

association_anova <- list()
association_fixed <- list()

for (band_name in levels(droplevels(combined$Band))) {
  for (measure in level_measures) {
    for (window_sec in levels(droplevels(combined$WindowSec))) {
      model_data <- combined[
        combined$Band == band_name &
          combined$WindowSec == window_sec &
          is.finite(combined[[measure]]) &
          is.finite(combined$z_Mean_LZs),
      ]
      model_data <- droplevels(model_data)
      model <- fit_lme_safely(model_data, paste0(measure, " ~ z_Mean_LZs + Eyes + Epoch"))
      extra <- list(
        WindowSec = as.character(window_sec),
        N = nrow(model_data),
        N_Participants = length(unique(model_data$Participant))
      )
      association_anova <- append_result(
        association_anova,
        anova_table(model, "PSD predicted by LZs + Eyes + Epoch", band_name, measure, extra)
      )
      association_fixed <- append_result(
        association_fixed,
        fixed_table(model, "PSD predicted by LZs + Eyes + Epoch", band_name, measure, extra)
      )
    }
  }
}

level_anova <- bind_results(level_anova)
level_fixed <- bind_results(level_fixed)
delta_anova <- bind_results(delta_anova)
delta_fixed <- bind_results(delta_fixed)
association_anova <- bind_results(association_anova)
association_fixed <- bind_results(association_fixed)

write.csv(level_anova, file.path(output_dir, "nlme_model_results_psd_level_anova.csv"), row.names = FALSE)
write.csv(level_fixed, file.path(output_dir, "nlme_model_results_psd_level_fixed_effects.csv"), row.names = FALSE)
write.csv(delta_anova, file.path(output_dir, "nlme_model_results_psd_delta_anova.csv"), row.names = FALSE)
write.csv(delta_fixed, file.path(output_dir, "nlme_model_results_psd_delta_fixed_effects.csv"), row.names = FALSE)
write.csv(association_anova, file.path(output_dir, "nlme_model_results_psd_lzs_association_anova.csv"), row.names = FALSE)
write.csv(association_fixed, file.path(output_dir, "nlme_model_results_psd_lzs_association_fixed_effects.csv"), row.names = FALSE)

cat("Wrote nlme mixed-model outputs:\n")
cat("  nlme_model_results_psd_level_anova.csv\n")
cat("  nlme_model_results_psd_level_fixed_effects.csv\n")
cat("  nlme_model_results_psd_delta_anova.csv\n")
cat("  nlme_model_results_psd_delta_fixed_effects.csv\n")
cat("  nlme_model_results_psd_lzs_association_anova.csv\n")
cat("  nlme_model_results_psd_lzs_association_fixed_effects.csv\n\n")

cat("Suggested first read:\n")
cat("  nlme_model_results_psd_level_anova.csv\n")
cat("  nlme_model_results_psd_delta_anova.csv\n")
cat("  nlme_model_results_psd_lzs_association_fixed_effects.csv\n")
