options(stringsAsFactors = FALSE)

base_dir <- "/Users/luisluna/Documents/Codex/2026-07-05/g/Mixed effects/outputs/PSD_band_LZs_mixed_models"
fig_dir <- file.path(base_dir, "figures", "poster")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

fixed_path <- file.path(base_dir, "alpha_psd_lzs_meq_fixed_effects.csv")
anova_path <- file.path(base_dir, "alpha_psd_lzs_meq_anova.csv")

fixed <- read.csv(fixed_path, check.names = FALSE)
anova_tbl <- read.csv(anova_path, check.names = FALSE)

poster_cex <- list(
  axis = 1.45,
  lab = 1.7,
  main = 1.85,
  sub = 1.25,
  annot = 1.35,
  legend = 1.35
)

sig_label <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < .001, "***",
      ifelse(p < .01, "**",
        ifelse(p < .05, "*", "ns"))))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < .001, "p < .001", sprintf("p = %.3f", p)))
}

beta_rows_for_window <- function(window_sec) {
  keep_lzc <- fixed$WindowSec == window_sec &
    fixed$ModelSet == "LZs_only" &
    fixed$ModelID %in% c("M1_simple", "M2_time_controlled") &
    fixed$Term == "z_Mean_LZs"
  keep_psd <- fixed$WindowSec == window_sec &
    fixed$ModelSet == "AlphaPSD_only" &
    fixed$ModelID %in% c("M1_simple", "M2_time_controlled") &
    fixed$Term == "z_Log10Power"
  keep_combined <- fixed$WindowSec == window_sec &
    fixed$ModelSet == "LZs_plus_AlphaPSD" &
    fixed$ModelID == "M2_time_controlled" &
    fixed$Term %in% c("z_Mean_LZs", "z_Log10Power")

  out <- fixed[keep_lzc | keep_psd | keep_combined, ]
  out$Predictor <- ifelse(out$Term == "z_Mean_LZs", "LZc", "Alpha PSD")
  out$ModelLabel <- ifelse(out$ModelID == "M1_simple", "Simple",
    ifelse(out$ModelSet == "LZs_plus_AlphaPSD", "combined\n+ time", "+ time"))
  out$Lower <- out$Estimate - 1.96 * out$Std_Error
  out$Upper <- out$Estimate + 1.96 * out$Std_Error
  out
}

interaction_rows_for_window <- function(window_sec) {
  keep <- anova_tbl$WindowSec == window_sec &
    anova_tbl$ModelID == "M3_time_interaction" &
    anova_tbl$Term %in% c("z_Mean_LZs:Epoch", "z_Log10Power:Epoch")

  out <- anova_tbl[keep, ]
  out$Predictor <- ifelse(out$Term == "z_Mean_LZs:Epoch", "LZc x time", "Alpha PSD x time")
  out$ModelLabel <- ifelse(out$ModelSet == "LZs_plus_AlphaPSD",
    ifelse(out$Predictor == "LZc x time", "Combined\nLZc x time", "Combined\nPSD x time"),
    out$Predictor)
  out$Evidence <- -log10(out$p_value)
  out
}

write_summary_csv <- function(window_sec, beta_tbl, interaction_tbl) {
  beta_out <- data.frame(
    WindowSec = window_sec,
    Panel = "Fixed-effect beta",
    Predictor = beta_tbl$Predictor,
    Model = beta_tbl$ModelLabel,
    Estimate = beta_tbl$Estimate,
    Lower95 = beta_tbl$Lower,
    Upper95 = beta_tbl$Upper,
    t = beta_tbl$t_value,
    df = beta_tbl$DF,
    p = beta_tbl$p_value,
    Significance = sig_label(beta_tbl$p_value)
  )

  int_out <- data.frame(
    WindowSec = window_sec,
    Panel = "Omnibus interaction with time",
    Predictor = interaction_tbl$Predictor,
    Model = interaction_tbl$ModelLabel,
    Estimate = interaction_tbl$F_value,
    Lower95 = NA,
    Upper95 = NA,
    t = NA,
    df = paste0(interaction_tbl$numDF, ",", interaction_tbl$denDF),
    p = interaction_tbl$p_value,
    Significance = sig_label(interaction_tbl$p_value)
  )

  write.csv(
    rbind(beta_out, int_out),
    file.path(base_dir, sprintf("poster_meq_all_four_model_summary_window%ss.csv", window_sec)),
    row.names = FALSE
  )
}

draw_errorbar <- function(x, y, low, high, col) {
  arrows(x, low, x, high, angle = 90, code = 3, length = 0.055, lwd = 3.5, col = col)
  points(x, y, pch = 19, cex = 1.7, col = col)
}

make_figure <- function(window_sec) {
  beta_tbl <- beta_rows_for_window(window_sec)
  interaction_tbl <- interaction_rows_for_window(window_sec)
  write_summary_csv(window_sec, beta_tbl, interaction_tbl)

  out_path <- file.path(fig_dir, sprintf("poster_meq_all_four_model_summary_window%ss.png", window_sec))
  png(out_path, width = 5000, height = 3300, res = 220)
  on.exit(dev.off(), add = TRUE)

  layout(matrix(c(1, 2), nrow = 1), widths = c(1.25, 1))
  par(family = "sans", fg = "black", col.axis = "black", col.lab = "black", col.main = "black")

  colors <- c("LZc" = "#d62728", "Alpha PSD" = "#1f77b4")

  par(mar = c(11.5, 11, 5.5, 2.5), xpd = NA)
  plot(NA, xlim = c(0.5, 6.5), ylim = c(-1.7, 1.8), xaxt = "n",
       xlab = "", ylab = "Fixed-effect estimate predicting MEQ4",
       main = "A. Additive predictor effects",
       cex.axis = poster_cex$axis, cex.lab = poster_cex$lab, cex.main = poster_cex$main)
  grid(col = "grey88", lty = 2)
  abline(h = 0, lty = 2, col = "grey45", lwd = 2)

  lzc <- beta_tbl[beta_tbl$Predictor == "LZc", ]
  psd <- beta_tbl[beta_tbl$Predictor == "Alpha PSD", ]
  x_lzc <- c(1, 2, 3)
  x_psd <- c(4.2, 5.2, 6.2)

  for (i in seq_len(nrow(lzc))) {
    draw_errorbar(x_lzc[i], lzc$Estimate[i], lzc$Lower[i], lzc$Upper[i], colors["LZc"])
    text(x_lzc[i], lzc$Upper[i] + .16, sig_label(lzc$p_value[i]), cex = poster_cex$annot, font = 2)
  }
  for (i in seq_len(nrow(psd))) {
    draw_errorbar(x_psd[i], psd$Estimate[i], psd$Lower[i], psd$Upper[i], colors["Alpha PSD"])
    text(x_psd[i], psd$Upper[i] + .16, sig_label(psd$p_value[i]), cex = poster_cex$annot, font = 2)
  }

  axis(1, at = c(x_lzc, x_psd),
       labels = c(lzc$ModelLabel, psd$ModelLabel),
       las = 2, cex.axis = poster_cex$axis)
  mtext("LZc predicting MEQ4", side = 3, line = 0.2, at = mean(x_lzc), col = colors["LZc"],
        cex = poster_cex$annot, font = 2)
  mtext("Alpha PSD predicting MEQ4", side = 3, line = 0.2, at = mean(x_psd), col = colors["Alpha PSD"],
        cex = poster_cex$annot, font = 2)
  par(mar = c(11.5, 10, 5.5, 2.5), xpd = NA)
  interaction_tbl <- interaction_tbl[order(interaction_tbl$ModelSet, interaction_tbl$Predictor), ]
  bar_cols <- ifelse(grepl("LZc", interaction_tbl$Predictor), colors["LZc"], colors["Alpha PSD"])
  bar_x <- seq_len(nrow(interaction_tbl))
  ylim_top <- max(1.8, interaction_tbl$Evidence, na.rm = TRUE) + 0.8
  plot(NA, xlim = c(0.4, length(bar_x) + 0.6), ylim = c(0, ylim_top), xaxt = "n",
       xlab = "", ylab = expression("Statistical evidence (" * -log[10](italic(p)) * ")"),
       main = "B. Omnibus predictor x time tests",
       cex.axis = poster_cex$axis, cex.lab = poster_cex$lab, cex.main = poster_cex$main)
  grid(col = "grey88", lty = 2)
  abline(h = -log10(.05), lty = 2, col = "grey45", lwd = 2)
  text(length(bar_x) + 0.25, -log10(.05) + 0.08, "p = .05", cex = poster_cex$sub, col = "grey35")
  rect(bar_x - .32, 0, bar_x + .32, interaction_tbl$Evidence, col = bar_cols, border = NA)
  text(bar_x, interaction_tbl$Evidence + .16,
       paste0(sig_label(interaction_tbl$p_value), "\n", fmt_p(interaction_tbl$p_value)),
       cex = poster_cex$sub, font = 2)
  axis(1, at = bar_x, labels = interaction_tbl$ModelLabel, las = 2, cex.axis = poster_cex$axis)

  invisible(out_path)
}

make_figure(2)
make_figure(10)

message("Created all-four MEQ4 model figures in: ", fig_dir)
