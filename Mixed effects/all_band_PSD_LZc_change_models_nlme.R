# =============================================================================
# OPTE all-band PSD change, LZc change, and MEQ4 models using base R + nlme.
#
# This extends the alpha-focused workflow to every PSD band:
#   delta, theta, alpha, beta, gamma
#
# Main question:
#   Does change in each PSD band from baseline line up with change in LZc from
#   baseline, after accounting for eye condition, dosing epoch, and repeated
#   measurements within participant?
#
# Secondary question:
#   Do PSD-band changes and LZc changes predict MEQ4 after accounting for
#   dosing epoch? MEQ4 has no eye condition, so those models use values averaged
#   across EC/EO within Participant x Epoch x Band x WindowSec.
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
fig_dir <- file.path(output_dir, "figures", "all_band_PSD_LZc_change")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

input_csv <- file.path(output_dir, "analysis_ready_psd_band_lzs_by_epoch_window.csv")
env_input_csv <- Sys.getenv("OPTE_PSD_LZS_READY_CSV", unset = "")
if (nzchar(env_input_csv)) input_csv <- env_input_csv

if (!file.exists(input_csv)) {
  stop(
    "Analysis-ready PSD + LZc CSV not found:\n",
    input_csv,
    "\nRun build_PSD_LZs_mixed_dataset_baseR.R first.",
    call. = FALSE
  )
}

cat("\nUsing PSD + LZc CSV:\n", input_csv, "\n", sep = "")
cat("Writing all-band PSD/LZc change outputs to:\n", output_dir, "\n\n", sep = "")

df <- read.csv(input_csv, stringsAsFactors = FALSE, colClasses = c(Participant = "character"))

required_cols <- c(
  "Participant", "Session", "Eyes", "Epoch", "Band", "WindowSec",
  "Log10Power", "RelativePower", "AbsolutePower",
  "Delta_Log10Power", "Delta_RelativePower", "Delta_AbsolutePower",
  "Mean_LZs", "Mean_LZsN", "z_Mean_LZs"
)
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Input table is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

clean_participant <- function(x) sprintf("%03d", as.integer(as.character(x)))
epoch_order_all <- paste0("Epoch", 0:4)
epoch_order_dosing <- paste0("Epoch", 1:4)
band_order <- c("delta", "theta", "alpha", "beta", "gamma")

df$Participant <- clean_participant(df$Participant)
df$Session <- tolower(trimws(as.character(df$Session)))
df$Eyes <- toupper(trimws(as.character(df$Eyes)))
df$Epoch <- trimws(as.character(df$Epoch))
df$Band <- tolower(trimws(as.character(df$Band)))
df$WindowSec <- as.integer(df$WindowSec)

# Compute LZc baseline-relative change. LZc is duplicated across bands after
# merging with PSD, so include Band in the merge to keep the table rectangular.
baseline_lzc <- df[df$Session == "baseline" & df$Epoch == "Epoch0",
                   c("Participant", "Eyes", "Band", "WindowSec", "Mean_LZs", "Mean_LZsN")]
names(baseline_lzc)[names(baseline_lzc) == "Mean_LZs"] <- "Baseline_Mean_LZs"
names(baseline_lzc)[names(baseline_lzc) == "Mean_LZsN"] <- "Baseline_Mean_LZsN"

df <- merge(
  df,
  baseline_lzc,
  by = c("Participant", "Eyes", "Band", "WindowSec"),
  all.x = TRUE,
  sort = FALSE
)
df$Delta_Mean_LZs <- df$Mean_LZs - df$Baseline_Mean_LZs
df$Delta_Mean_LZsN <- df$Mean_LZsN - df$Baseline_Mean_LZsN

z_by_group <- function(data, value_col, group_cols, out_col) {
  data[[out_col]] <- NA_real_
  key <- interaction(data[group_cols], drop = TRUE, lex.order = TRUE)
  for (level in levels(key)) {
    idx <- which(key == level)
    values <- data[[value_col]][idx]
    if (sum(is.finite(values)) >= 2 && is.finite(sd(values, na.rm = TRUE)) && sd(values, na.rm = TRUE) > 0) {
      data[[out_col]][idx] <- as.numeric(scale(values))
    }
  }
  data
}

df <- z_by_group(df, "Delta_Mean_LZs", c("Participant", "WindowSec"), "z_Delta_Mean_LZs")
df <- z_by_group(df, "Delta_Mean_LZsN", c("Participant", "WindowSec"), "z_Delta_Mean_LZsN")

write.csv(df, file.path(output_dir, "all_band_psd_lzc_by_eye_window.csv"), row.names = FALSE)

mean_aggregate <- function(data, value_cols, key_cols) {
  out <- NULL
  for (value_col in value_cols) {
    form <- as.formula(paste(value_col, "~", paste(key_cols, collapse = " + ")))
    tmp <- aggregate(form, data = data, FUN = function(x) mean(x, na.rm = TRUE))
    if (is.null(out)) {
      out <- tmp
    } else {
      out <- merge(out, tmp, by = key_cols, all = TRUE, sort = FALSE)
    }
  }
  out
}

epoch_avg <- mean_aggregate(
  df,
  c(
    "Log10Power", "RelativePower", "AbsolutePower",
    "Delta_Log10Power", "Delta_RelativePower", "Delta_AbsolutePower",
    "Mean_LZs", "Mean_LZsN", "Delta_Mean_LZs", "Delta_Mean_LZsN",
    "z_Mean_LZs", "z_Delta_Mean_LZs"
  ),
  c("Participant", "Session", "Epoch", "Band", "WindowSec")
)
epoch_avg <- epoch_avg[order(epoch_avg$Participant, epoch_avg$Band, epoch_avg$Epoch, epoch_avg$WindowSec), ]
write.csv(epoch_avg, file.path(output_dir, "all_band_psd_lzc_epoch_avg_for_meq.csv"), row.names = FALSE)

meq_data <- data.frame(
  Participant = clean_participant(c(5, 7, 8, 11, 12, 13, 16, 17, 18, 19)),
  T4 = c(1, 3, NA, 0, 5, 3, 3, 5, 2, 5),
  T6 = c(3, 3, 5, 1, 5, 5, 3, 5, 4, 3),
  T8 = c(1, 4, 2, 1, 3, 0, 1, 5, 4, 0),
  T9 = c(1, 2, 2, 0, 3, 0, 1, 4, 4, 0),
  stringsAsFactors = FALSE
)
meq_long <- rbind(
  data.frame(Participant = meq_data$Participant, MEQ_timepoint = "T4", Epoch = "Epoch1", MEQ4 = meq_data$T4),
  data.frame(Participant = meq_data$Participant, MEQ_timepoint = "T6", Epoch = "Epoch2", MEQ4 = meq_data$T6),
  data.frame(Participant = meq_data$Participant, MEQ_timepoint = "T8", Epoch = "Epoch3", MEQ4 = meq_data$T8),
  data.frame(Participant = meq_data$Participant, MEQ_timepoint = "T9", Epoch = "Epoch4", MEQ4 = meq_data$T9)
)

meq_epoch <- merge(
  epoch_avg[epoch_avg$Session == "dosing" & epoch_avg$Epoch != "Epoch0", ],
  meq_long,
  by = c("Participant", "Epoch"),
  all.x = TRUE,
  sort = FALSE
)
meq_epoch <- meq_epoch[is.finite(meq_epoch$MEQ4), ]
write.csv(meq_epoch, file.path(output_dir, "all_band_psd_lzc_meq_epoch_avg.csv"), row.names = FALSE)

correlation_row <- function(data, x_col, y_col, label, band, window_sec = NA, eyes = NA, dosing_only = FALSE) {
  keep <- is.finite(data[[x_col]]) & is.finite(data[[y_col]])
  tmp <- data[keep, ]
  if (nrow(tmp) < 4) {
    return(data.frame(
      Analysis = label, Band = band, X = x_col, Y = y_col, WindowSec = window_sec, Eyes = eyes,
      DosingOnly = dosing_only, N = nrow(tmp), r = NA_real_, p_value = NA_real_,
      CI95_low = NA_real_, CI95_high = NA_real_, stringsAsFactors = FALSE
    ))
  }
  ct <- cor.test(tmp[[x_col]], tmp[[y_col]], method = "pearson")
  ci <- suppressWarnings(ct$conf.int)
  data.frame(
    Analysis = label, Band = band, X = x_col, Y = y_col, WindowSec = window_sec, Eyes = eyes,
    DosingOnly = dosing_only, N = nrow(tmp), r = unname(ct$estimate), p_value = ct$p.value,
    CI95_low = ci[1], CI95_high = ci[2], stringsAsFactors = FALSE
  )
}

fit_lme_safely <- function(data, formula_text, epoch_levels) {
  data$Participant <- factor(data$Participant)
  if ("Eyes" %in% names(data)) data$Eyes <- factor(data$Eyes, levels = c("EC", "EO"))
  data$Epoch <- factor(data$Epoch, levels = epoch_levels)
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

fixed_effect_rows <- function(model, model_name, band, window_sec, outcome, predictors, n, n_participants) {
  if (inherits(model, "error")) {
    return(data.frame(
      Model = model_name, Band = band, WindowSec = window_sec, Outcome = outcome, Predictors = predictors,
      Term = NA_character_, Estimate = NA_real_, Std_Error = NA_real_, DF = NA_real_,
      t_value = NA_real_, p_value = NA_real_, N = n, N_Participants = n_participants,
      Error = model$message, stringsAsFactors = FALSE
    ))
  }
  coefs <- as.data.frame(summary(model)$tTable)
  coefs$Term <- rownames(coefs)
  rownames(coefs) <- NULL
  data.frame(
    Model = model_name, Band = band, WindowSec = window_sec, Outcome = outcome, Predictors = predictors,
    Term = coefs$Term, Estimate = coefs$Value, Std_Error = coefs$Std.Error, DF = coefs$DF,
    t_value = coefs$`t-value`, p_value = coefs$`p-value`, N = n, N_Participants = n_participants,
    Error = NA_character_, stringsAsFactors = FALSE
  )
}

cor_rows <- list()
model_rows <- list()
meq_model_rows <- list()

primary_delta_measures <- c("Delta_Log10Power", "Delta_RelativePower", "Delta_AbsolutePower")

for (band in band_order) {
  for (window_sec in sort(unique(df$WindowSec))) {
    for (eye in sort(unique(df$Eyes))) {
      sub <- df[df$Band == band & df$WindowSec == window_sec & df$Eyes == eye, ]
      dosing_sub <- sub[sub$Session == "dosing" & sub$Epoch != "Epoch0", ]
      for (measure in primary_delta_measures) {
        cor_rows[[length(cor_rows) + 1]] <- correlation_row(
          dosing_sub, measure, "Delta_Mean_LZs",
          "PSD-band change vs LZc change by eye",
          band, window_sec, eye, TRUE
        )
      }
    }

    avg_sub <- epoch_avg[epoch_avg$Band == band & epoch_avg$WindowSec == window_sec, ]
    avg_dosing_sub <- avg_sub[avg_sub$Session == "dosing" & avg_sub$Epoch != "Epoch0", ]
    for (measure in primary_delta_measures) {
      cor_rows[[length(cor_rows) + 1]] <- correlation_row(
        avg_dosing_sub, measure, "Delta_Mean_LZs",
        "PSD-band change vs LZc change averaged across eyes",
        band, window_sec, "EC+EO average", TRUE
      )
    }

    for (measure in primary_delta_measures) {
      model_data <- df[
        df$Band == band & df$WindowSec == window_sec &
          df$Session == "dosing" & df$Epoch != "Epoch0" &
          is.finite(df[[measure]]) & is.finite(df$Delta_Mean_LZs),
      ]
      model <- fit_lme_safely(
        model_data,
        paste0(measure, " ~ Delta_Mean_LZs + Eyes + Epoch"),
        epoch_order_dosing
      )
      model_rows[[length(model_rows) + 1]] <- fixed_effect_rows(
        model,
        "PSD-band change predicted by LZc change + Eyes + Epoch",
        band, window_sec, measure, "Delta_Mean_LZs + Eyes + Epoch",
        nrow(model_data), length(unique(model_data$Participant))
      )
    }

    meq_data_band <- meq_epoch[
      meq_epoch$Band == band & meq_epoch$WindowSec == window_sec &
        is.finite(meq_epoch$MEQ4) & is.finite(meq_epoch$Delta_Log10Power) & is.finite(meq_epoch$Delta_Mean_LZs),
    ]
    meq_model <- fit_lme_safely(
      meq_data_band,
      "MEQ4 ~ Delta_Log10Power + Delta_Mean_LZs + Epoch",
      epoch_order_dosing
    )
    meq_model_rows[[length(meq_model_rows) + 1]] <- fixed_effect_rows(
      meq_model,
      "MEQ4 predicted by PSD-band change + LZc change + Epoch",
      band, window_sec, "MEQ4", "Delta_Log10Power + Delta_Mean_LZs + Epoch",
      nrow(meq_data_band), length(unique(meq_data_band$Participant))
    )
  }
}

correlations <- do.call(rbind, cor_rows)
models <- do.call(rbind, model_rows)
meq_models <- do.call(rbind, meq_model_rows)

write.csv(correlations, file.path(output_dir, "all_band_psd_lzc_change_correlations.csv"), row.names = FALSE)
write.csv(models, file.path(output_dir, "all_band_psd_lzc_change_nlme_fixed_effects.csv"), row.names = FALSE)
write.csv(meq_models, file.path(output_dir, "all_band_meq_psd_lzc_change_nlme_fixed_effects.csv"), row.names = FALSE)

format_p <- function(p) {
  if (!is.finite(p)) return("p = NA")
  if (p < 0.001) return("p < .001")
  paste0("p = ", sub("^0", "", sprintf("%.3f", p)))
}

band_cols <- c(delta = "#1f77b4", theta = "#ff7f0e", alpha = "#d62728", beta = "#9467bd", gamma = "#2ca02c")

for (window_sec in sort(unique(epoch_avg$WindowSec))) {
  png(file.path(fig_dir, paste0("all_band_psd_lzc_change_scatter_window", window_sec, "s.png")),
      width = 1900, height = 1200, res = 180, bg = "white")
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(bg = "white", fg = "black", col.axis = "black", col.lab = "black", col.main = "black",
      mfrow = c(2, 3), mar = c(4.6, 4.8, 3.0, 1.2), oma = c(0, 0, 2.0, 0))

  for (band in band_order) {
    plot_df <- epoch_avg[
      epoch_avg$Band == band & epoch_avg$WindowSec == window_sec &
        epoch_avg$Session == "dosing" & epoch_avg$Epoch != "Epoch0" &
        is.finite(epoch_avg$Delta_Mean_LZs) & is.finite(epoch_avg$Delta_Log10Power),
    ]
    if (nrow(plot_df) < 4) {
      plot.new()
      title(main = band)
      text(0.5, 0.5, "No data")
      next
    }
    plot(
      plot_df$Delta_Mean_LZs, plot_df$Delta_Log10Power,
      pch = 19, cex = 1.1, col = band_cols[band],
      xlab = "LZc change from baseline",
      ylab = "Log10 power change from baseline",
      main = paste0(toupper(substr(band, 1, 1)), substr(band, 2, nchar(band)))
    )
    grid(col = "gray88")
    abline(h = 0, v = 0, lty = 2, col = "gray55")
    abline(lm(Delta_Log10Power ~ Delta_Mean_LZs, data = plot_df), col = "black", lwd = 2)
    cor_row <- correlations[
      correlations$Band == band &
        correlations$WindowSec == window_sec &
        correlations$Eyes == "EC+EO average" &
        correlations$X == "Delta_Log10Power",
    ]
    if (nrow(cor_row) > 0) {
      legend(
        "topright",
        legend = sprintf("r = %.2f, %s", cor_row$r[1], format_p(cor_row$p_value[1])),
        bty = "n", cex = 1.05
      )
    }
  }
  plot.new()
  title(main = "Notes")
  text(0.05, 0.75, "Each panel averages EC/EO within participant and epoch.", adj = 0, cex = 1.105)
  text(0.05, 0.55, "Mixed-effects tables include eye condition and epoch controls.", adj = 0, cex = 1.105)
  title(main = paste0("PSD-band change vs LZc change, window ", window_sec, "s"), outer = TRUE)
}

cat("Wrote all-band PSD/LZc change outputs:\n")
cat("  all_band_psd_lzc_by_eye_window.csv\n")
cat("  all_band_psd_lzc_epoch_avg_for_meq.csv\n")
cat("  all_band_psd_lzc_meq_epoch_avg.csv\n")
cat("  all_band_psd_lzc_change_correlations.csv\n")
cat("  all_band_psd_lzc_change_nlme_fixed_effects.csv\n")
cat("  all_band_meq_psd_lzc_change_nlme_fixed_effects.csv\n")
cat("Figures written to:\n", fig_dir, "\n", sep = "")
