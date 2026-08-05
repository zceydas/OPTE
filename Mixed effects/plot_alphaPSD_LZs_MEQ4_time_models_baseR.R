# =============================================================================
# Figures for MEQ4 time models with LZs and alpha PSD.
#
# Run after:
#   Rscript "Mixed effects/LZs_alphaPSD_predicts_MEQ4_time_model_nlme.R"
#
# Figures emphasize the same question as the original reference analysis script:
#   - Do EEG metrics predict MEQ4 before time is included?
#   - What happens after Epoch/time is added?
#   - Do LZs and alpha PSD add independent information above time?
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

dataset_csv <- file.path(output_dir, "alpha_psd_lzs_meq_model_dataset.csv")
fixed_csv <- file.path(output_dir, "alpha_psd_lzs_meq_fixed_effects.csv")
increment_csv <- file.path(output_dir, "alpha_psd_lzs_meq_incremental_lrt.csv")

for (path in c(dataset_csv, fixed_csv, increment_csv)) {
  if (!file.exists(path)) {
    stop(
      "Required input not found:\n",
      path,
      "\nRun LZs_alphaPSD_predicts_MEQ4_time_model_nlme.R first.",
      call. = FALSE
    )
  }
}

df <- read.csv(dataset_csv, stringsAsFactors = FALSE, colClasses = c(Participant = "character"))
fixed <- read.csv(fixed_csv, stringsAsFactors = FALSE)
increment <- read.csv(increment_csv, stringsAsFactors = FALSE)

epoch_levels <- paste0("Epoch", 1:4)
epoch_labels <- c(Epoch1 = "Epoch1\nT4", Epoch2 = "Epoch2\nT6", Epoch3 = "Epoch3\nT8", Epoch4 = "Epoch4\nT9")
epoch_cols <- c(Epoch1 = "#1f77b4", Epoch2 = "#2ca02c", Epoch3 = "#d62728", Epoch4 = "#9467bd")

theme_par <- function() {
  par(
    bg = "white", fg = "black", col.axis = "black", col.lab = "black",
    col.main = "black", family = "sans", las = 1, bty = "l", cex.lab = 1.70, cex.axis = 1.45, cex.main = 1.55
  )
}

safe_dev_off <- function() {
  if (names(dev.cur()) != "null device") dev.off()
}

mean_ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) return(c(N = 0, Mean = NA_real_, SE = NA_real_, CI95_low = NA_real_, CI95_high = NA_real_))
  se <- if (n > 1) sd(x) / sqrt(n) else NA_real_
  ci <- if (n > 1) qt(0.975, df = n - 1) * se else NA_real_
  c(N = n, Mean = mean(x), SE = se, CI95_low = mean(x) - ci, CI95_high = mean(x) + ci)
}

sig_label <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

metric_label <- function(term) {
  ifelse(term == "z_Mean_LZs", "LZs", "Alpha PSD")
}

model_label <- function(model_id, model_set) {
  if (model_id == "M1_simple") return("Simple")
  if (model_id == "M2_time_controlled" && model_set %in% c("LZs_only", "AlphaPSD_only")) return("+ Epoch")
  if (model_id == "M2_time_controlled" && model_set == "LZs_plus_AlphaPSD") return("Combined\n+ Epoch")
  if (model_id == "M3_time_interaction") return("Interaction\nwith Epoch")
  model_id
}

plot_scatter_panels <- function(sub, window_sec) {
  png(
    file.path(fig_dir, paste0("meq4_simple_scatter_window", window_sec, "s.png")),
    width = 1800, height = 900, res = 180, bg = "white"
  )
  on.exit(safe_dev_off(), add = TRUE)
  theme_par()
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

  panels <- list(
    list(x = "z_Mean_LZs", xlab = "LZs z-score", main = "MEQ4 vs LZs"),
    list(x = "z_Log10Power", xlab = "Alpha log10 PSD z-score", main = "MEQ4 vs alpha PSD")
  )
  for (panel in panels) {
    x <- sub[[panel$x]]
    y <- sub$MEQ4
    keep <- is.finite(x) & is.finite(y)
    cols <- epoch_cols[sub$Epoch[keep]]
    plot(
      x[keep], y[keep], pch = 19, cex = 1.25, col = cols,
      xlab = panel$xlab, ylab = "MEQ4", main = paste0(panel$main, " (", window_sec, "s window)")
    )
    grid(col = "gray88")
    abline(lm(y[keep] ~ x[keep]), lwd = 2, col = "black")
    ct <- cor.test(x[keep], y[keep], method = "pearson")
    usr <- par("usr")
    text(
      usr[1] + 0.04 * diff(usr[1:2]), usr[4] - 0.08 * diff(usr[3:4]),
      labels = sprintf("r = %.2f, p = %.3f", unname(ct$estimate), ct$p.value),
      adj = c(0, 1), cex = 1.20
    )
    legend(
      "bottomright", legend = names(epoch_cols), col = epoch_cols, pch = 19,
      bty = "n", cex = 1.10
    )
  }
}

plot_coefficient_attenuation <- function(fixed_sub, window_sec) {
  keep <- fixed_sub$Term %in% c("z_Mean_LZs", "z_Log10Power") &
    fixed_sub$ModelID %in% c("M1_simple", "M2_time_controlled") &
    fixed_sub$ModelSet %in% c("LZs_only", "AlphaPSD_only", "LZs_plus_AlphaPSD")
  coefs <- fixed_sub[keep, ]
  coefs <- coefs[!(coefs$ModelSet == "LZs_plus_AlphaPSD" & coefs$ModelID == "M1_simple"), ]
  coefs$CI95_low <- coefs$Estimate - 1.96 * coefs$Std_Error
  coefs$CI95_high <- coefs$Estimate + 1.96 * coefs$Std_Error
  coefs$Label <- mapply(model_label, coefs$ModelID, coefs$ModelSet)
  coefs$Metric <- metric_label(coefs$Term)

  order_keys <- c(
    "LZs|Simple", "LZs|+ Epoch", "LZs|Combined\n+ Epoch",
    "Alpha PSD|Simple", "Alpha PSD|+ Epoch", "Alpha PSD|Combined\n+ Epoch"
  )
  coefs$Key <- paste(coefs$Metric, coefs$Label, sep = "|")
  coefs <- coefs[match(order_keys, coefs$Key, nomatch = 0), ]

  png(
    file.path(fig_dir, paste0("meq4_predictor_attenuation_window", window_sec, "s.png")),
    width = 1800, height = 1000, res = 180, bg = "white"
  )
  on.exit(safe_dev_off(), add = TRUE)
  theme_par()
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(7, 5, 4, 2))

  x <- seq_len(nrow(coefs))
  cols <- ifelse(coefs$Metric == "LZs", "#d62728", "#1f77b4")
  ylim <- range(c(coefs$CI95_low, coefs$CI95_high, 0), na.rm = TRUE)
  ylim <- ylim + c(-0.12, 0.20) * diff(ylim)
  plot(
    x, coefs$Estimate, type = "n", xaxt = "n", ylim = ylim,
    xlab = "", ylab = "Fixed-effect estimate predicting MEQ4",
    main = paste0("MEQ4 predictor effects before and after time control (", window_sec, "s LZs window)")
  )
  grid(col = "gray88")
  abline(h = 0, lty = 2, col = "gray45")
  segments(x, coefs$CI95_low, x, coefs$CI95_high, col = cols, lwd = 2)
  points(x, coefs$Estimate, pch = 19, cex = 1.4, col = cols)
  axis(1, at = x, labels = paste(coefs$Metric, coefs$Label, sep = "\n"), las = 2, cex.axis = 1.05)
  text(
    x, coefs$CI95_high + 0.04 * diff(ylim),
    labels = vapply(coefs$p_value, sig_label, character(1)),
    cex = 1.20
  )
  legend(
    "topright", legend = c("LZs", "Alpha PSD"),
    col = c("#d62728", "#1f77b4"), pch = 19, bty = "n"
  )
  mtext("*** p < .001, ** p < .01, * p < .05, ns p >= .05", side = 1, line = 5.8, cex = 1.05)
}

plot_epoch_trajectories <- function(sub, window_sec) {
  sub$z_MEQ4 <- as.numeric(scale(sub$MEQ4))
  vars <- c("z_Mean_LZs", "z_Log10Power", "z_MEQ4")
  labels <- c(z_Mean_LZs = "LZs", z_Log10Power = "Alpha PSD", z_MEQ4 = "MEQ4")
  cols <- c(z_Mean_LZs = "#d62728", z_Log10Power = "#1f77b4", z_MEQ4 = "#2ca02c")

  rows <- list()
  for (epoch in epoch_levels) {
    epoch_sub <- sub[sub$Epoch == epoch, ]
    for (var in vars) {
      stats <- mean_ci(epoch_sub[[var]])
      rows[[length(rows) + 1]] <- data.frame(
        Epoch = epoch, Variable = var, N = stats["N"], Mean = stats["Mean"],
        CI95_low = stats["CI95_low"], CI95_high = stats["CI95_high"],
        stringsAsFactors = FALSE
      )
    }
  }
  summary_df <- do.call(rbind, rows)

  png(
    file.path(fig_dir, paste0("meq4_epoch_trajectories_window", window_sec, "s.png")),
    width = 1800, height = 1000, res = 180, bg = "white"
  )
  on.exit(safe_dev_off(), add = TRUE)
  theme_par()
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(5, 5, 4, 2))

  ylim <- range(c(summary_df$CI95_low, summary_df$CI95_high), na.rm = TRUE)
  plot(
    seq_along(epoch_levels), rep(NA_real_, length(epoch_levels)),
    xaxt = "n", ylim = ylim, xlab = "", ylab = "Standardized mean +/- 95% CI",
    main = paste0("Dosing-epoch trajectories used in time-controlled MEQ4 models (", window_sec, "s)")
  )
  grid(col = "gray88")
  axis(1, at = seq_along(epoch_levels), labels = epoch_labels[epoch_levels])

  for (var in vars) {
    tmp <- summary_df[summary_df$Variable == var, ]
    x <- match(tmp$Epoch, epoch_levels)
    lines(x, tmp$Mean, col = cols[var], lwd = 2)
    points(x, tmp$Mean, col = cols[var], pch = 19, cex = 1.2)
    segments(x, tmp$CI95_low, x, tmp$CI95_high, col = cols[var], lwd = 1.6)
  }
  legend(
    "topright", legend = labels[vars], col = cols[vars],
    lwd = 2, pch = 19, bty = "n"
  )
}

plot_incremental_lrt <- function(increment_sub, window_sec) {
  full_rows <- increment_sub[!is.na(increment_sub$Test) & increment_sub$Test != "", ]
  if (nrow(full_rows) == 0) return(invisible(FALSE))
  labels <- ifelse(
    grepl("alpha PSD add", full_rows$ModelSet),
    "Alpha PSD added\nto LZs + Epoch",
    "LZs added\nto alpha PSD + Epoch"
  )
  pvals <- full_rows$p.value
  lrt <- full_rows$L.Ratio

  png(
    file.path(fig_dir, paste0("meq4_incremental_lrt_window", window_sec, "s.png")),
    width = 1500, height = 1000, res = 180, bg = "white"
  )
  on.exit(safe_dev_off(), add = TRUE)
  theme_par()
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(6, 5, 4, 2))

  cols <- c("#1f77b4", "#d62728")
  ylim <- c(0, max(lrt, na.rm = TRUE) * 1.35)
  bp <- barplot(
    lrt, names.arg = labels, col = cols, border = NA,
    ylim = ylim,
    ylab = "Likelihood-ratio statistic",
    main = paste0("Incremental MEQ4 model tests after controlling Epoch (", window_sec, "s)")
  )
  abline(h = 0, col = "gray45")
  grid(nx = NA, ny = NULL, col = "gray88")
  text(
    bp, lrt + 0.08 * max(lrt, na.rm = TRUE),
    labels = sprintf("p = %.3f\n%s", pvals, vapply(pvals, sig_label, character(1))),
    cex = 1.20
  )
}

for (window_sec in sort(unique(df$WindowSec))) {
  sub <- df[df$WindowSec == window_sec, ]
  fixed_sub <- fixed[fixed$WindowSec == window_sec, ]
  increment_sub <- increment[increment$WindowSec == window_sec, ]

  plot_scatter_panels(sub, window_sec)
  plot_coefficient_attenuation(fixed_sub, window_sec)
  plot_epoch_trajectories(sub, window_sec)
  plot_incremental_lrt(increment_sub, window_sec)
}

cat("Wrote MEQ4 figures to:\n", fig_dir, "\n\n", sep = "")
cat("Figure files:\n")
cat("  meq4_simple_scatter_window2s.png\n")
cat("  meq4_simple_scatter_window10s.png\n")
cat("  meq4_predictor_attenuation_window2s.png\n")
cat("  meq4_predictor_attenuation_window10s.png\n")
cat("  meq4_epoch_trajectories_window2s.png\n")
cat("  meq4_epoch_trajectories_window10s.png\n")
cat("  meq4_incremental_lrt_window2s.png\n")
cat("  meq4_incremental_lrt_window10s.png\n")
