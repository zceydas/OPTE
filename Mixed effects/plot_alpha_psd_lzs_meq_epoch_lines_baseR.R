# =============================================================================
# Alpha PSD, LZs, and MEQ4 epoch line graphs with error bars.
#
# Goal:
#   Make the main results easy to inspect visually across dosing epochs:
#   - mean alpha PSD, LZs, and MEQ4 trajectories
#   - epoch-by-epoch correlations between alpha PSD, LZs, and MEQ4
#
# Inputs are produced by:
#   alpha_PSD_LZs_MEQ_models_nlme.R
#
# Outputs are PNG line graphs in:
#   Mixed effects/outputs/PSD_band_LZs_mixed_models/figures
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
alpha_meq_csv <- file.path(output_dir, "alpha_psd_lzs_meq_epoch_avg.csv")

if (!file.exists(alpha_eye_csv) || !file.exists(alpha_meq_csv)) {
  stop(
    "Alpha PSD-LZs-MEQ CSVs not found. Run alpha_PSD_LZs_MEQ_models_nlme.R first.",
    call. = FALSE
  )
}

cat("\nUsing alpha eye-level CSV:\n", alpha_eye_csv, "\n", sep = "")
cat("Using alpha MEQ CSV:\n", alpha_meq_csv, "\n", sep = "")
cat("Writing epoch line figures to:\n", fig_dir, "\n\n", sep = "")

alpha_eye <- read.csv(alpha_eye_csv, stringsAsFactors = FALSE, colClasses = c(Participant = "character"))
alpha_meq <- read.csv(alpha_meq_csv, stringsAsFactors = FALSE, colClasses = c(Participant = "character"))

epoch_order <- paste0("Epoch", 1:4)
epoch_labels <- c("Epoch1\nT4", "Epoch2\nT6", "Epoch3\nT8", "Epoch4\nT9")
epoch_order_with_baseline <- paste0("Epoch", 0:4)
epoch_labels_with_baseline <- c("Baseline\nEpoch0", epoch_labels)
window_order <- sort(unique(alpha_meq$WindowSec))

z_all <- function(x) {
  if (sum(is.finite(x)) < 2 || sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
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

cor_ci <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  if (n < 4 || sd(x) == 0 || sd(y) == 0) {
    return(c(N = n, r = NA_real_, CI95_low = NA_real_, CI95_high = NA_real_, p_value = NA_real_))
  }
  ct <- cor.test(x, y, method = "pearson")
  c(N = n, r = unname(ct$estimate), CI95_low = ct$conf.int[1], CI95_high = ct$conf.int[2], p_value = ct$p.value)
}

sig_label <- function(p) {
  vapply(p, function(one_p) {
    if (!is.finite(one_p)) return("")
    if (one_p < 0.001) return("***")
    if (one_p < 0.01) return("**")
    if (one_p < 0.05) return("*")
    "ns"
  }, character(1))
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

line_cols <- c(
  "Alpha PSD" = "#1f77b4",
  "LZs" = "#d62728",
  "MEQ4" = "#2ca02c",
  "PSD-LZs" = "#7f3c8d",
  "PSD-MEQ4" = "#11a579",
  "LZs-MEQ4" = "#3969ac",
  "EC" = "#1f77b4",
  "EO" = "#d62728",
  "EC+EO average" = "#2ca02c"
)

plot_epoch_lines <- function(summary_df, title_text, ylab_text, file_name, local_epoch_order = epoch_order, local_epoch_labels = epoch_labels) {
  y_range <- range(c(summary_df$Mean - 1.96 * summary_df$SE, summary_df$Mean + 1.96 * summary_df$SE), na.rm = TRUE)
  y_pad <- max(0.05, diff(y_range) * 0.08)
  x_values <- seq_along(local_epoch_order)
  write_plot(file_name)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mar = c(5.5, 5.2, 4.2, 2))
  plot(
    NA,
    xlim = c(0.8, length(local_epoch_order) + 0.2),
    ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
    xaxt = "n",
    xlab = "",
    ylab = ylab_text,
    main = title_text,
    cex.lab = 1.70,
    cex.axis = 1.45,
    cex.main = 1.55
  )
  axis(1, at = x_values, labels = local_epoch_labels, cex.axis = 1.20)
  grid(nx = NA, ny = NULL, col = "gray88")
  abline(h = 0, col = "gray55", lty = 2)
  if ("Epoch0" %in% local_epoch_order) {
    abline(v = 1.5, col = "gray70", lty = 2, lwd = 1.1)
  }

  for (signal in unique(summary_df$Signal)) {
    sub <- summary_df[summary_df$Signal == signal, ]
    y <- rep(NA_real_, length(local_epoch_order))
    se <- rep(NA_real_, length(local_epoch_order))
    idx <- match(sub$Epoch, local_epoch_order)
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
          col = line_cols[[signal]], lwd = 1.2
        )
      }
      lines(x_values, y, type = "b", pch = 19, lwd = 2.4, cex = 1, col = line_cols[[signal]])
    }
  }
  legend("topright", legend = unique(summary_df$Signal), col = line_cols[unique(summary_df$Signal)], lty = 1, pch = 19, lwd = 2.4, bty = "n")
}

plot_correlation_lines <- function(cor_df, title_text, file_name, local_epoch_order = epoch_order, local_epoch_labels = epoch_labels) {
  y_range <- range(c(cor_df$CI95_low, cor_df$CI95_high), na.rm = TRUE)
  y_pad <- max(0.05, diff(y_range) * 0.08)
  x_values <- seq_along(local_epoch_order)
  write_plot(file_name)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mar = c(5.5, 5.2, 4.2, 2))
  plot(
    NA,
    xlim = c(0.8, length(local_epoch_order) + 0.2),
    ylim = c(max(-1, y_range[1] - y_pad), min(1, y_range[2] + y_pad)),
    xaxt = "n",
    xlab = "",
    ylab = "Pearson r",
    main = title_text,
    cex.lab = 1.70,
    cex.axis = 1.45,
    cex.main = 1.55
  )
  axis(1, at = x_values, labels = local_epoch_labels, cex.axis = 1.20)
  grid(nx = NA, ny = NULL, col = "gray88")
  abline(h = 0, col = "gray55", lty = 2)
  if ("Epoch0" %in% local_epoch_order) {
    abline(v = 1.5, col = "gray70", lty = 2, lwd = 1.1)
  }

  pair_order <- unique(cor_df$Pair)
  pair_offsets <- setNames(seq(-0.08, 0.08, length.out = length(pair_order)), pair_order)

  for (pair in pair_order) {
    sub <- cor_df[cor_df$Pair == pair, ]
    y <- rep(NA_real_, length(local_epoch_order))
    lo <- rep(NA_real_, length(local_epoch_order))
    hi <- rep(NA_real_, length(local_epoch_order))
    idx <- match(sub$Epoch, local_epoch_order)
    y[idx] <- sub$r
    lo[idx] <- sub$CI95_low
    hi[idx] <- sub$CI95_high
    finite <- is.finite(y)
    if (any(finite)) {
      x_pos <- x_values + pair_offsets[[pair]]
      arrows(
        x_pos[finite], lo[finite],
        x_pos[finite], hi[finite],
        angle = 90, code = 3, length = 0.035,
        col = line_cols[[pair]], lwd = 1.2
      )
      lines(x_pos, y, type = "b", pch = 19, lwd = 2.4, cex = 1, col = line_cols[[pair]])
      matched <- match(local_epoch_order[finite], sub$Epoch)
      labels <- paste0(
        sig_label(sub$p_value[matched]),
        "\n",
        "n=", sub$N[matched]
      )
      text(x_pos[finite], y[finite], labels = labels, pos = 3, cex = 1.15, col = "black")
    }
  }
  legend("bottomleft", legend = unique(cor_df$Pair), col = line_cols[unique(cor_df$Pair)], lty = 1, pch = 19, lwd = 2.4, bty = "n")
  legend("topright", legend = c("*** p < .001", "** p < .01", "* p < .05", "ns p >= .05"), bty = "n", cex = 1.05)
}

all_summary_rows <- list()
all_cor_rows <- list()
all_baseline_summary_rows <- list()
all_baseline_cor_rows <- list()

for (window_sec in window_order) {
  baseline_sub <- alpha_eye[
    alpha_eye$WindowSec == window_sec &
      alpha_eye$Epoch %in% epoch_order_with_baseline,
  ]
  baseline_sub$z_Log10Power_all <- z_all(baseline_sub$Log10Power)
  baseline_sub$z_Mean_LZs_all <- z_all(baseline_sub$Mean_LZs)

  baseline_long <- rbind(
    data.frame(Participant = baseline_sub$Participant, Epoch = baseline_sub$Epoch, WindowSec = window_sec, Signal = "Alpha PSD", Value = baseline_sub$z_Log10Power_all),
    data.frame(Participant = baseline_sub$Participant, Epoch = baseline_sub$Epoch, WindowSec = window_sec, Signal = "LZs", Value = baseline_sub$z_Mean_LZs_all)
  )
  baseline_summary <- mean_se(baseline_long, "Value", c("WindowSec", "Signal", "Epoch"))
  all_baseline_summary_rows[[length(all_baseline_summary_rows) + 1]] <- baseline_summary
  plot_epoch_lines(
    baseline_summary,
    paste0("Alpha PSD and LZs From Baseline Through Dosing Epochs (", window_sec, "s LZs window)"),
    "Standardized value (mean +/- 95% CI)",
    paste0("alpha_psd_lzs_with_baseline_epoch_lines_window", window_sec, "s.png"),
    epoch_order_with_baseline,
    epoch_labels_with_baseline
  )

  baseline_cor_rows <- list()
  for (epoch in epoch_order_with_baseline) {
    epoch_df <- baseline_sub[baseline_sub$Epoch == epoch, ]
    a <- cor_ci(epoch_df$z_Log10Power_all, epoch_df$z_Mean_LZs_all)
    baseline_cor_rows[[length(baseline_cor_rows) + 1]] <- data.frame(WindowSec = window_sec, Epoch = epoch, Pair = "PSD-LZs", t(a), stringsAsFactors = FALSE)
  }
  baseline_cor_summary <- do.call(rbind, baseline_cor_rows)
  numeric_cols <- c("N", "r", "CI95_low", "CI95_high", "p_value")
  for (col in numeric_cols) baseline_cor_summary[[col]] <- as.numeric(baseline_cor_summary[[col]])
  all_baseline_cor_rows[[length(all_baseline_cor_rows) + 1]] <- baseline_cor_summary
  plot_correlation_lines(
    baseline_cor_summary,
    paste0("Baseline-to-Dosing Correlations: Alpha PSD and LZs (", window_sec, "s LZs window)"),
    paste0("alpha_psd_lzs_with_baseline_epoch_correlations_window", window_sec, "s.png"),
    epoch_order_with_baseline,
    epoch_labels_with_baseline
  )

  alpha_epoch_all <- aggregate(
    cbind(Log10Power, Mean_LZs, Delta_Log10Power, Delta_Mean_LZs) ~ Participant + Session + Epoch + WindowSec,
    data = baseline_sub,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  alpha_epoch_all <- alpha_epoch_all[alpha_epoch_all$Epoch %in% epoch_order_with_baseline, ]
  meq_for_window <- unique(alpha_meq[alpha_meq$WindowSec == window_sec, c("Participant", "Epoch", "MEQ4")])
  sub <- merge(alpha_epoch_all, meq_for_window, by = c("Participant", "Epoch"), all.x = TRUE, sort = FALSE)
  sub <- sub[order(sub$Participant, match(sub$Epoch, epoch_order_with_baseline)), ]
  sub$z_MEQ4 <- z_all(sub$MEQ4)
  sub$Signal_Alpha_PSD <- z_all(sub$Log10Power)
  sub$Signal_LZs <- z_all(sub$Mean_LZs)
  sub$Signal_MEQ4 <- sub$z_MEQ4

  long_signal <- rbind(
    data.frame(Participant = sub$Participant, Epoch = sub$Epoch, WindowSec = window_sec, Signal = "Alpha PSD", Value = sub$Signal_Alpha_PSD),
    data.frame(Participant = sub$Participant, Epoch = sub$Epoch, WindowSec = window_sec, Signal = "LZs", Value = sub$Signal_LZs),
    data.frame(Participant = sub$Participant, Epoch = sub$Epoch, WindowSec = window_sec, Signal = "MEQ4", Value = sub$Signal_MEQ4)
  )
  signal_summary <- mean_se(long_signal, "Value", c("WindowSec", "Signal", "Epoch"))
  names(signal_summary)[names(signal_summary) == "Value"] <- "Value"
  all_summary_rows[[length(all_summary_rows) + 1]] <- signal_summary
  plot_epoch_lines(
    signal_summary,
    paste0("Alpha PSD, LZs, and MEQ4 From Baseline Through Dosing Epochs (", window_sec, "s LZs window)"),
    "Standardized value (mean +/- 95% CI)",
    paste0("alpha_psd_lzs_meq_epoch_lines_window", window_sec, "s.png"),
    epoch_order_with_baseline,
    epoch_labels_with_baseline
  )

  sub_delta <- sub
  sub_delta$z_MEQ4 <- z_all(sub_delta$MEQ4)
  sub_delta$Signal_Alpha_PSD <- z_all(sub_delta$Delta_Log10Power)
  sub_delta$Signal_LZs <- z_all(sub_delta$Delta_Mean_LZs)
  sub_delta$Signal_MEQ4 <- sub_delta$z_MEQ4
  long_delta <- rbind(
    data.frame(Participant = sub_delta$Participant, Epoch = sub_delta$Epoch, WindowSec = window_sec, Signal = "Alpha PSD", Value = sub_delta$Signal_Alpha_PSD),
    data.frame(Participant = sub_delta$Participant, Epoch = sub_delta$Epoch, WindowSec = window_sec, Signal = "LZs", Value = sub_delta$Signal_LZs),
    data.frame(Participant = sub_delta$Participant, Epoch = sub_delta$Epoch, WindowSec = window_sec, Signal = "MEQ4", Value = sub_delta$Signal_MEQ4)
  )
  delta_summary <- mean_se(long_delta, "Value", c("WindowSec", "Signal", "Epoch"))
  all_summary_rows[[length(all_summary_rows) + 1]] <- transform(delta_summary, Signal = paste0(Signal, " change"))
  plot_epoch_lines(
    delta_summary,
    paste0("Alpha PSD Change, LZs Change, and MEQ4 From Baseline Through Dosing Epochs (", window_sec, "s LZs window)"),
    "Standardized value (mean +/- 95% CI)",
    paste0("alpha_psd_lzs_meq_delta_epoch_lines_window", window_sec, "s.png"),
    epoch_order_with_baseline,
    epoch_labels_with_baseline
  )

  cor_rows <- list()
  for (epoch in epoch_order_with_baseline) {
    epoch_df <- sub[sub$Epoch == epoch, ]
    a <- cor_ci(epoch_df$Signal_Alpha_PSD, epoch_df$Signal_LZs)
    cor_rows[[length(cor_rows) + 1]] <- data.frame(WindowSec = window_sec, Epoch = epoch, Pair = "PSD-LZs", t(a), stringsAsFactors = FALSE)
    a <- cor_ci(epoch_df$Signal_Alpha_PSD, epoch_df$MEQ4)
    cor_rows[[length(cor_rows) + 1]] <- data.frame(WindowSec = window_sec, Epoch = epoch, Pair = "PSD-MEQ4", t(a), stringsAsFactors = FALSE)
    a <- cor_ci(epoch_df$Signal_LZs, epoch_df$MEQ4)
    cor_rows[[length(cor_rows) + 1]] <- data.frame(WindowSec = window_sec, Epoch = epoch, Pair = "LZs-MEQ4", t(a), stringsAsFactors = FALSE)
  }
  cor_summary <- do.call(rbind, cor_rows)
  numeric_cols <- c("N", "r", "CI95_low", "CI95_high", "p_value")
  for (col in numeric_cols) cor_summary[[col]] <- as.numeric(cor_summary[[col]])
  all_cor_rows[[length(all_cor_rows) + 1]] <- cor_summary
  plot_correlation_lines(
    cor_summary,
    paste0("Baseline-to-Dosing Correlations: Alpha PSD, LZs, and MEQ4 (", window_sec, "s LZs window)"),
    paste0("alpha_psd_lzs_meq_epoch_correlations_window", window_sec, "s.png"),
    epoch_order_with_baseline,
    epoch_labels_with_baseline
  )
}

trajectory_summary <- do.call(rbind, all_summary_rows)
correlation_summary <- do.call(rbind, all_cor_rows)
baseline_trajectory_summary <- do.call(rbind, all_baseline_summary_rows)
baseline_correlation_summary <- do.call(rbind, all_baseline_cor_rows)
write.csv(trajectory_summary, file.path(output_dir, "alpha_psd_lzs_meq_epoch_line_summary.csv"), row.names = FALSE)
write.csv(correlation_summary, file.path(output_dir, "alpha_psd_lzs_meq_epoch_correlations.csv"), row.names = FALSE)
write.csv(baseline_trajectory_summary, file.path(output_dir, "alpha_psd_lzs_with_baseline_epoch_line_summary.csv"), row.names = FALSE)
write.csv(baseline_correlation_summary, file.path(output_dir, "alpha_psd_lzs_with_baseline_epoch_correlations.csv"), row.names = FALSE)

cat("Wrote epoch line summary:\n")
cat("  alpha_psd_lzs_meq_epoch_line_summary.csv\n")
cat("  alpha_psd_lzs_meq_epoch_correlations.csv\n\n")
cat("  alpha_psd_lzs_with_baseline_epoch_line_summary.csv\n")
cat("  alpha_psd_lzs_with_baseline_epoch_correlations.csv\n\n")
cat("Wrote line figures:\n")
for (window_sec in window_order) {
  cat("  ", paste0("alpha_psd_lzs_with_baseline_epoch_lines_window", window_sec, "s.png"), "\n", sep = "")
  cat("  ", paste0("alpha_psd_lzs_with_baseline_epoch_correlations_window", window_sec, "s.png"), "\n", sep = "")
  cat("  ", paste0("alpha_psd_lzs_meq_epoch_lines_window", window_sec, "s.png"), "\n", sep = "")
  cat("  ", paste0("alpha_psd_lzs_meq_delta_epoch_lines_window", window_sec, "s.png"), "\n", sep = "")
  cat("  ", paste0("alpha_psd_lzs_meq_epoch_correlations_window", window_sec, "s.png"), "\n", sep = "")
}
