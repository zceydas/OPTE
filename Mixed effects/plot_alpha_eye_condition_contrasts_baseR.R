# =============================================================================
# Alpha EEG eye-condition line graphs with error bars.
#
# Outputs:
#   1. Alpha PSD: EC vs EO, baseline through dosing
#   2. LZs: EC vs EO, baseline through dosing
#   3. Alpha PSD: EO minus EC contrast, baseline through dosing
#   4. LZs: EO minus EC contrast, baseline through dosing
#
# Inputs are produced by:
#   alpha_PSD_LZs_MEQ_models_nlme.R
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
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

alpha_eye_csv <- file.path(output_dir, "alpha_psd_lzs_by_eye_window.csv")
if (!file.exists(alpha_eye_csv)) {
  stop(
    "Alpha eye-level CSV not found. Run alpha_PSD_LZs_MEQ_models_nlme.R first.",
    call. = FALSE
  )
}

cat("\nUsing alpha eye-level CSV:\n", alpha_eye_csv, "\n", sep = "")
cat("Writing eye-condition figures to:\n", fig_dir, "\n\n", sep = "")

alpha_eye <- read.csv(alpha_eye_csv, stringsAsFactors = FALSE, colClasses = c(Participant = "character"))

required_cols <- c("Participant", "Eyes", "WindowSec", "Session", "Epoch", "Log10Power", "Mean_LZs")
missing_cols <- setdiff(required_cols, names(alpha_eye))
if (length(missing_cols) > 0) {
  stop("Alpha eye-level table is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

alpha_eye$Participant <- sprintf("%03d", as.integer(as.character(alpha_eye$Participant)))
alpha_eye$Eyes <- toupper(trimws(as.character(alpha_eye$Eyes)))
alpha_eye$Epoch <- trimws(as.character(alpha_eye$Epoch))
alpha_eye$WindowSec <- as.integer(alpha_eye$WindowSec)

epoch_order <- paste0("Epoch", 0:4)
epoch_labels <- c("Baseline\nEpoch0", "Epoch1\nT4", "Epoch2\nT6", "Epoch3\nT8", "Epoch4\nT9")
window_order <- sort(unique(alpha_eye$WindowSec))

line_cols <- c(
  "EC" = "#1f77b4",
  "EO" = "#d62728",
  "EO-EC" = "#7f3c8d"
)

sig_label <- function(p) {
  vapply(p, function(one_p) {
    if (!is.finite(one_p)) return("")
    if (one_p < 0.001) return("***")
    if (one_p < 0.01) return("**")
    if (one_p < 0.05) return("*")
    "ns"
  }, character(1))
}

mean_se <- function(df, value_col, group_cols) {
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
  out$key_tmp <- NULL
  out
}

paired_eye_p <- function(df, value_col) {
  rows <- list()
  for (epoch in epoch_order) {
    sub <- df[df$Epoch == epoch, ]
    wide <- reshape(
      sub[, c("Participant", "Eyes", value_col)],
      idvar = "Participant",
      timevar = "Eyes",
      direction = "wide"
    )
    ec_col <- paste0(value_col, ".EC")
    eo_col <- paste0(value_col, ".EO")
    keep <- is.finite(wide[[ec_col]]) & is.finite(wide[[eo_col]])
    p <- NA_real_
    n <- sum(keep)
    if (n >= 3) {
      p <- t.test(wide[[eo_col]][keep], wide[[ec_col]][keep], paired = TRUE)$p.value
    }
    rows[[length(rows) + 1]] <- data.frame(Epoch = epoch, N = n, p_value = p, Sig = sig_label(p), stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

one_sample_p <- function(df, value_col) {
  rows <- list()
  for (epoch in epoch_order) {
    sub <- df[df$Epoch == epoch, ]
    values <- sub[[value_col]]
    values <- values[is.finite(values)]
    p <- NA_real_
    n <- length(values)
    if (n >= 3 && sd(values) > 0) {
      p <- t.test(values, mu = 0)$p.value
    }
    rows[[length(rows) + 1]] <- data.frame(Epoch = epoch, N = n, p_value = p, Sig = sig_label(p), stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

write_plot <- function(file_name, width = 1800, height = 1200) {
  png(file.path(fig_dir, file_name), width = width, height = height, res = 180, bg = "white")
  par(
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
}

plot_eye_lines <- function(summary_df, stats_df, title_text, ylab_text, file_name) {
  y_range <- range(c(summary_df$Mean - 1.96 * summary_df$SE, summary_df$Mean + 1.96 * summary_df$SE), na.rm = TRUE)
  y_pad <- max(0.05, diff(y_range) * 0.12)
  x_values <- seq_along(epoch_order)
  write_plot(file_name)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mar = c(5.5, 5.2, 4.2, 2))
  plot(
    NA,
    xlim = c(0.8, length(epoch_order) + 0.2),
    ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
    xaxt = "n",
    xlab = "",
    ylab = ylab_text,
    main = title_text,
    cex.lab = 1.70,
    cex.axis = 1.45,
    cex.main = 1.55
  )
  axis(1, at = x_values, labels = epoch_labels, cex.axis = 1.20)
  grid(nx = NA, ny = NULL, col = "gray88")
  abline(v = 1.5, col = "gray70", lty = 2, lwd = 1.1)

  for (eye in c("EC", "EO")) {
    sub <- summary_df[summary_df$Eyes == eye, ]
    y <- rep(NA_real_, length(epoch_order))
    se <- rep(NA_real_, length(epoch_order))
    idx <- match(sub$Epoch, epoch_order)
    y[idx] <- sub$Mean
    se[idx] <- sub$SE
    finite <- is.finite(y)
    if (any(finite)) {
      err_idx <- finite & is.finite(se) & se > 0
      if (any(err_idx)) {
        arrows(
          x_values[err_idx], y[err_idx] - 1.96 * se[err_idx],
          x_values[err_idx], y[err_idx] + 1.96 * se[err_idx],
          angle = 90, code = 3, length = 0.035,
          col = line_cols[[eye]], lwd = 1.2
        )
      }
      lines(x_values, y, type = "b", pch = 19, lwd = 2.4, cex = 1, col = line_cols[[eye]])
    }
  }

  for (i in seq_len(nrow(stats_df))) {
    x <- match(stats_df$Epoch[i], epoch_order)
    if (is.finite(x) && nzchar(stats_df$Sig[i])) {
      text(x, y_range[2] + y_pad * 0.55, labels = paste0(stats_df$Sig[i], "\nn=", stats_df$N[i]), cex = 1.15, col = "black")
    }
  }

  legend("topright", legend = c("EC", "EO"), col = line_cols[c("EC", "EO")], lty = 1, pch = 19, lwd = 2.4, bty = "n")
  legend("bottomleft", legend = c("*** p < .001", "** p < .01", "* p < .05", "ns p >= .05"), bty = "n", cex = 1.05)
}

plot_contrast_line <- function(summary_df, stats_df, title_text, ylab_text, file_name) {
  y_range <- range(c(summary_df$Mean - 1.96 * summary_df$SE, summary_df$Mean + 1.96 * summary_df$SE), na.rm = TRUE)
  y_pad <- max(0.02, diff(y_range) * 0.12)
  x_values <- seq_along(epoch_order)
  write_plot(file_name)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mar = c(5.5, 5.2, 4.2, 2))
  plot(
    NA,
    xlim = c(0.8, length(epoch_order) + 0.2),
    ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
    xaxt = "n",
    xlab = "",
    ylab = ylab_text,
    main = title_text,
    cex.lab = 1.70,
    cex.axis = 1.45,
    cex.main = 1.55
  )
  axis(1, at = x_values, labels = epoch_labels, cex.axis = 1.20)
  grid(nx = NA, ny = NULL, col = "gray88")
  abline(h = 0, col = "gray55", lty = 2)
  abline(v = 1.5, col = "gray70", lty = 2, lwd = 1.1)

  y <- rep(NA_real_, length(epoch_order))
  se <- rep(NA_real_, length(epoch_order))
  idx <- match(summary_df$Epoch, epoch_order)
  y[idx] <- summary_df$Mean
  se[idx] <- summary_df$SE
  finite <- is.finite(y)
  if (any(finite)) {
    err_idx <- finite & is.finite(se) & se > 0
    if (any(err_idx)) {
      arrows(
        x_values[err_idx], y[err_idx] - 1.96 * se[err_idx],
        x_values[err_idx], y[err_idx] + 1.96 * se[err_idx],
        angle = 90, code = 3, length = 0.035,
        col = line_cols[["EO-EC"]], lwd = 1.2
      )
    }
    lines(x_values, y, type = "b", pch = 19, lwd = 2.4, cex = 1, col = line_cols[["EO-EC"]])
  }

  for (i in seq_len(nrow(stats_df))) {
    x <- match(stats_df$Epoch[i], epoch_order)
    if (is.finite(x) && nzchar(stats_df$Sig[i])) {
      text(x, y_range[2] + y_pad * 0.55, labels = paste0(stats_df$Sig[i], "\nn=", stats_df$N[i]), cex = 1.15, col = "black")
    }
  }

  legend("topright", legend = "EO - EC", col = line_cols[["EO-EC"]], lty = 1, pch = 19, lwd = 2.4, bty = "n")
  legend("bottomleft", legend = c("*** p < .001", "** p < .01", "* p < .05", "ns p >= .05"), bty = "n", cex = 1.05)
}

all_eye_summaries <- list()
all_eye_stats <- list()
all_contrast_summaries <- list()
all_contrast_stats <- list()

for (window_sec in window_order) {
  sub <- alpha_eye[alpha_eye$WindowSec == window_sec & alpha_eye$Epoch %in% epoch_order, ]

  psd_summary <- mean_se(sub, "Log10Power", c("WindowSec", "Eyes", "Epoch"))
  psd_summary$Measure <- "Alpha_Log10Power"
  psd_stats <- paired_eye_p(sub, "Log10Power")
  psd_stats$WindowSec <- window_sec
  psd_stats$Measure <- "Alpha_Log10Power_EO_vs_EC"
  all_eye_summaries[[length(all_eye_summaries) + 1]] <- psd_summary
  all_eye_stats[[length(all_eye_stats) + 1]] <- psd_stats
  plot_eye_lines(
    psd_summary,
    psd_stats,
    paste0("Alpha PSD by Eye Condition: Baseline Through Dosing (", window_sec, "s LZs window)"),
    "Alpha log10 power (mean +/- 95% CI)",
    paste0("alpha_psd_eye_condition_lines_window", window_sec, "s.png")
  )

  lzs_summary <- mean_se(sub, "Mean_LZs", c("WindowSec", "Eyes", "Epoch"))
  lzs_summary$Measure <- "Mean_LZs"
  lzs_stats <- paired_eye_p(sub, "Mean_LZs")
  lzs_stats$WindowSec <- window_sec
  lzs_stats$Measure <- "Mean_LZs_EO_vs_EC"
  all_eye_summaries[[length(all_eye_summaries) + 1]] <- lzs_summary
  all_eye_stats[[length(all_eye_stats) + 1]] <- lzs_stats
  plot_eye_lines(
    lzs_summary,
    lzs_stats,
    paste0("LZs by Eye Condition: Baseline Through Dosing (", window_sec, "s LZs window)"),
    "Mean LZs (mean +/- 95% CI)",
    paste0("lzs_eye_condition_lines_window", window_sec, "s.png")
  )

  wide_psd <- reshape(sub[, c("Participant", "Epoch", "Eyes", "Log10Power", "Mean_LZs")], idvar = c("Participant", "Epoch"), timevar = "Eyes", direction = "wide")
  wide_psd$Delta_Log10Power_EO_minus_EC <- wide_psd$Log10Power.EO - wide_psd$Log10Power.EC
  wide_psd$Delta_Mean_LZs_EO_minus_EC <- wide_psd$Mean_LZs.EO - wide_psd$Mean_LZs.EC
  wide_psd$WindowSec <- window_sec

  psd_contrast_summary <- mean_se(wide_psd, "Delta_Log10Power_EO_minus_EC", c("WindowSec", "Epoch"))
  psd_contrast_summary$Measure <- "Delta_Log10Power_EO_minus_EC"
  psd_contrast_stats <- one_sample_p(wide_psd, "Delta_Log10Power_EO_minus_EC")
  psd_contrast_stats$WindowSec <- window_sec
  psd_contrast_stats$Measure <- "Delta_Log10Power_EO_minus_EC"
  all_contrast_summaries[[length(all_contrast_summaries) + 1]] <- psd_contrast_summary
  all_contrast_stats[[length(all_contrast_stats) + 1]] <- psd_contrast_stats
  plot_contrast_line(
    psd_contrast_summary,
    psd_contrast_stats,
    paste0("Alpha PSD Eye Contrast: EO Minus EC (", window_sec, "s LZs window)"),
    "EO - EC alpha log10 power",
    paste0("alpha_psd_eo_minus_ec_lines_window", window_sec, "s.png")
  )

  lzs_contrast_summary <- mean_se(wide_psd, "Delta_Mean_LZs_EO_minus_EC", c("WindowSec", "Epoch"))
  lzs_contrast_summary$Measure <- "Delta_Mean_LZs_EO_minus_EC"
  lzs_contrast_stats <- one_sample_p(wide_psd, "Delta_Mean_LZs_EO_minus_EC")
  lzs_contrast_stats$WindowSec <- window_sec
  lzs_contrast_stats$Measure <- "Delta_Mean_LZs_EO_minus_EC"
  all_contrast_summaries[[length(all_contrast_summaries) + 1]] <- lzs_contrast_summary
  all_contrast_stats[[length(all_contrast_stats) + 1]] <- lzs_contrast_stats
  plot_contrast_line(
    lzs_contrast_summary,
    lzs_contrast_stats,
    paste0("LZs Eye Contrast: EO Minus EC (", window_sec, "s LZs window)"),
    "EO - EC Mean LZs",
    paste0("lzs_eo_minus_ec_lines_window", window_sec, "s.png")
  )
}

eye_summary <- do.call(rbind, all_eye_summaries)
eye_stats <- do.call(rbind, all_eye_stats)
contrast_summary <- do.call(rbind, all_contrast_summaries)
contrast_stats <- do.call(rbind, all_contrast_stats)

write.csv(eye_summary, file.path(output_dir, "alpha_eye_condition_line_summary.csv"), row.names = FALSE)
write.csv(eye_stats, file.path(output_dir, "alpha_eye_condition_paired_stats.csv"), row.names = FALSE)
write.csv(contrast_summary, file.path(output_dir, "alpha_eo_minus_ec_line_summary.csv"), row.names = FALSE)
write.csv(contrast_stats, file.path(output_dir, "alpha_eo_minus_ec_stats.csv"), row.names = FALSE)

cat("Wrote eye-condition summaries:\n")
cat("  alpha_eye_condition_line_summary.csv\n")
cat("  alpha_eye_condition_paired_stats.csv\n")
cat("  alpha_eo_minus_ec_line_summary.csv\n")
cat("  alpha_eo_minus_ec_stats.csv\n\n")

cat("Wrote eye-condition figures:\n")
for (window_sec in window_order) {
  cat("  ", paste0("alpha_psd_eye_condition_lines_window", window_sec, "s.png"), "\n", sep = "")
  cat("  ", paste0("lzs_eye_condition_lines_window", window_sec, "s.png"), "\n", sep = "")
  cat("  ", paste0("alpha_psd_eo_minus_ec_lines_window", window_sec, "s.png"), "\n", sep = "")
  cat("  ", paste0("lzs_eo_minus_ec_lines_window", window_sec, "s.png"), "\n", sep = "")
}
