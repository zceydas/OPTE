# =============================================================================
# Descriptive OPTE PSD band summaries using base R only.
#
# Purpose:
#   Summarize how delta/theta/alpha/beta/gamma power changes across baseline,
#   dosing epochs, and eye condition before running mixed-effects models.
#
# Inputs:
#   Defaults to the analysis-ready PSD table created by:
#   build_PSD_LZs_mixed_dataset_baseR.R
#
# Outputs:
#   - descriptive_psd_band_summary_by_epoch.csv
#   - descriptive_psd_band_delta_summary_by_epoch.csv
#   - PNG figures for absolute, log10, relative, and baseline-delta power
# =============================================================================

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
fig_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

psd_ready_csv <- file.path(output_dir, "analysis_ready_psd_band_by_epoch.csv")
env_psd_ready_csv <- Sys.getenv("OPTE_PSD_READY_CSV", unset = "")
if (nzchar(env_psd_ready_csv)) {
  psd_ready_csv <- env_psd_ready_csv
}

if (!file.exists(psd_ready_csv)) {
  stop(
    "Analysis-ready PSD CSV not found:\n",
    psd_ready_csv,
    "\nRun build_PSD_LZs_mixed_dataset_baseR.R first.",
    call. = FALSE
  )
}

cat("\nUsing analysis-ready PSD CSV:\n", psd_ready_csv, "\n", sep = "")
cat("Writing descriptive outputs to:\n", output_dir, "\n\n", sep = "")

psd <- read.csv(psd_ready_csv, stringsAsFactors = FALSE, colClasses = c(Participant = "character"))

required_cols <- c(
  "Participant", "Session", "Eyes", "Epoch", "Band",
  "AbsolutePower", "Log10Power", "RelativePower",
  "Delta_AbsolutePower", "Delta_Log10Power", "Delta_RelativePower"
)
missing_cols <- setdiff(required_cols, names(psd))
if (length(missing_cols) > 0) {
  stop("PSD table is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

band_order <- c("delta", "theta", "alpha", "beta", "gamma")
epoch_order <- paste0("Epoch", 0:4)
epoch_labels <- c("Baseline\nEpoch0", "Dosing\nEpoch1", "Dosing\nEpoch2", "Dosing\nEpoch3", "Dosing\nEpoch4")
eyes_order <- c("EC", "EO")
eyes_labels <- c(EC = "Eyes Closed", EO = "Eyes Open")
band_colors <- c(
  delta = "#1f77b4",
  theta = "#2ca02c",
  alpha = "#d62728",
  beta = "#9467bd",
  gamma = "#ff7f0e"
)

psd$Band <- factor(tolower(psd$Band), levels = band_order)
psd$Epoch <- factor(psd$Epoch, levels = epoch_order)
psd$Eyes <- factor(toupper(psd$Eyes), levels = eyes_order)
psd$Session <- tolower(psd$Session)

mean_se_table <- function(df, value_col, group_cols) {
  key <- interaction(df[group_cols], drop = TRUE, lex.order = TRUE)
  out <- unique(df[group_cols])
  out$key_tmp <- interaction(out[group_cols], drop = TRUE, lex.order = TRUE)
  values_by_key <- split(df[[value_col]], key)
  participants_by_key <- split(df$Participant, key)
  out$N <- as.integer(vapply(values_by_key[as.character(out$key_tmp)], function(x) sum(is.finite(x)), numeric(1)))
  out$N_Participants <- as.integer(vapply(participants_by_key[as.character(out$key_tmp)], function(x) length(unique(x)), numeric(1)))
  out$Mean <- as.numeric(vapply(values_by_key[as.character(out$key_tmp)], function(x) mean(x, na.rm = TRUE), numeric(1)))
  out$SD <- as.numeric(vapply(values_by_key[as.character(out$key_tmp)], function(x) sd(x, na.rm = TRUE), numeric(1)))
  out$SE <- out$SD / sqrt(out$N)
  out$CI95_low <- out$Mean - 1.96 * out$SE
  out$CI95_high <- out$Mean + 1.96 * out$SE
  out$Measure <- value_col
  out$key_tmp <- NULL
  out[order(out$Band, out$Eyes, out$Epoch), c(group_cols, "Measure", "N", "N_Participants", "Mean", "SD", "SE", "CI95_low", "CI95_high")]
}

level_measures <- c("AbsolutePower", "Log10Power", "RelativePower")
delta_measures <- c("Delta_AbsolutePower", "Delta_Log10Power", "Delta_RelativePower")

level_summary <- do.call(
  rbind,
  lapply(level_measures, function(measure) mean_se_table(psd, measure, c("Band", "Eyes", "Epoch")))
)

delta_summary <- do.call(
  rbind,
  lapply(delta_measures, function(measure) mean_se_table(psd, measure, c("Band", "Eyes", "Epoch")))
)

write.csv(level_summary, file.path(output_dir, "descriptive_psd_band_summary_by_epoch.csv"), row.names = FALSE)
write.csv(delta_summary, file.path(output_dir, "descriptive_psd_band_delta_summary_by_epoch.csv"), row.names = FALSE)

format_axis_label <- function(measure) {
  switch(
    measure,
    AbsolutePower = "Absolute power",
    Log10Power = "Log10 power",
    RelativePower = "Relative power",
    Delta_AbsolutePower = "Change from baseline: absolute power",
    Delta_Log10Power = "Change from baseline: log10 power",
    Delta_RelativePower = "Change from baseline: relative power",
    measure
  )
}

plot_measure <- function(summary_df, measure, file_name) {
  plot_df <- summary_df[summary_df$Measure == measure, ]
  y_range <- range(c(plot_df$CI95_low, plot_df$CI95_high), na.rm = TRUE)
  if (!all(is.finite(y_range)) || diff(y_range) == 0) {
    y_range <- range(plot_df$Mean, na.rm = TRUE)
  }
  y_pad <- max(1e-6, diff(y_range) * 0.08)
  y_lim <- c(y_range[1] - y_pad, y_range[2] + y_pad)

  png(file.path(fig_dir, file_name), width = 2200, height = 1200, res = 180, bg = "white")
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  par(
    mfrow = c(1, 2),
    mar = c(7.5, 5.2, 4.5, 1.5),
    oma = c(0, 0, 3.5, 0),
    bg = "white",
    fg = "black",
    col.axis = "black",
    col.lab = "black",
    col.main = "black",
    family = "sans",
    cex.lab = 1.70,
    cex.axis = 1.45,
    cex.main = 1.55
  )

  for (eye in eyes_order) {
    eye_df <- plot_df[plot_df$Eyes == eye, ]
    plot(
      NA,
      xlim = c(0.8, length(epoch_order) + 0.2),
      ylim = y_lim,
      xaxt = "n",
      xlab = "",
      ylab = format_axis_label(measure),
      main = eyes_labels[[eye]],
      cex.lab = 1.70,
      cex.axis = 1.45,
      cex.main = 1.55
    )
    axis(1, at = seq_along(epoch_order), labels = epoch_labels, las = 1, cex.axis = 1.15)
    grid(nx = NA, ny = NULL, col = "gray88", lty = 1)
    abline(v = 1.5, col = "gray70", lty = 2, lwd = 1.2)
    if (grepl("^Delta_", measure)) {
      abline(h = 0, col = "gray45", lty = 2, lwd = 1.2)
    }

    for (band in band_order) {
      band_df <- eye_df[eye_df$Band == band, ]
      y <- rep(NA_real_, length(epoch_order))
      se <- rep(NA_real_, length(epoch_order))
      idx <- match(as.character(band_df$Epoch), epoch_order)
      y[idx] <- band_df$Mean
      se[idx] <- band_df$SE
      x <- seq_along(epoch_order)
      finite <- is.finite(y)
      if (any(finite)) {
        errorbar_idx <- finite & is.finite(se) & se > 0
        if (any(errorbar_idx)) {
          arrows(
            x[errorbar_idx], y[errorbar_idx] - 1.96 * se[errorbar_idx],
            x[errorbar_idx], y[errorbar_idx] + 1.96 * se[errorbar_idx],
            angle = 90, code = 3, length = 0.035,
            col = band_colors[[band]], lwd = 1
          )
        }
        lines(x[finite], y[finite], type = "b", pch = 19, lwd = 2.2, cex = 1.20, col = band_colors[[band]])
      }
    }

    legend(
      "topright",
      legend = band_order,
      col = band_colors[band_order],
      lty = 1,
      pch = 19,
      lwd = 2,
      bty = "n",
      cex = 1.20
    )
  }

  title(main = paste("PSD Band Changes Across Baseline, Dosing Epochs, and Eye Condition:", format_axis_label(measure)), outer = TRUE, cex.main = 1.55)
}

plot_measure(level_summary, "AbsolutePower", "PSD_band_absolute_power_by_epoch_and_eyes.png")
plot_measure(level_summary, "Log10Power", "PSD_band_log10_power_by_epoch_and_eyes.png")
plot_measure(level_summary, "RelativePower", "PSD_band_relative_power_by_epoch_and_eyes.png")
plot_measure(delta_summary, "Delta_AbsolutePower", "PSD_band_delta_absolute_power_by_epoch_and_eyes.png")
plot_measure(delta_summary, "Delta_Log10Power", "PSD_band_delta_log10_power_by_epoch_and_eyes.png")
plot_measure(delta_summary, "Delta_RelativePower", "PSD_band_delta_relative_power_by_epoch_and_eyes.png")

cat("Wrote descriptive CSVs:\n")
cat("  descriptive_psd_band_summary_by_epoch.csv\n")
cat("  descriptive_psd_band_delta_summary_by_epoch.csv\n\n")
cat("Wrote figures to:\n", fig_dir, "\n\n", sep = "")
cat("Best first figures:\n")
cat("  PSD_band_log10_power_by_epoch_and_eyes.png\n")
cat("  PSD_band_relative_power_by_epoch_and_eyes.png\n")
cat("  PSD_band_delta_log10_power_by_epoch_and_eyes.png\n")
cat("  PSD_band_delta_relative_power_by_epoch_and_eyes.png\n")
