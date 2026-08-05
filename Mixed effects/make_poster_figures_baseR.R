# =============================================================================
# Poster-quality figures using the newer analysis logic.
#
# This script is intentionally written in readable sections rather than compact
# loops. It replaces the earlier draft figure scripts for poster purposes.
#
# Design choices:
#   - White background, black text.
#   - Legends outside the plotting region whenever possible.
#   - Large, high-contrast significance labels.
#   - No labels placed directly on top of data points.
#
# Inputs:
#   outputs/PSD_band_LZs_mixed_models/alpha_psd_lzs_by_eye_window.csv
#   outputs/PSD_band_LZs_mixed_models/alpha_psd_lzs_meq_epoch_avg.csv
#   outputs/PSD_band_LZs_mixed_models/alpha_psd_lzs_nlme_fixed_effects.csv
#   outputs/PSD_band_LZs_mixed_models/alpha_psd_lzs_meq_fixed_effects.csv
#   outputs/PSD_band_LZs_mixed_models/alpha_psd_lzs_meq_incremental_lrt.csv
#
# Outputs:
#   outputs/PSD_band_LZs_mixed_models/figures/poster/
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
fig_dir <- file.path(output_dir, "figures", "poster")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

alpha_eye_csv <- file.path(output_dir, "alpha_psd_lzs_by_eye_window.csv")
alpha_meq_csv <- file.path(output_dir, "alpha_psd_lzs_meq_epoch_avg.csv")
alpha_fixed_csv <- file.path(output_dir, "alpha_psd_lzs_nlme_fixed_effects.csv")
fixed_csv <- file.path(output_dir, "alpha_psd_lzs_meq_fixed_effects.csv")
increment_csv <- file.path(output_dir, "alpha_psd_lzs_meq_incremental_lrt.csv")

for (path in c(alpha_eye_csv, alpha_meq_csv, alpha_fixed_csv, fixed_csv, increment_csv)) {
  if (!file.exists(path)) {
    stop(
      "Required input not found:\n",
      path,
      "\nRun alpha_PSD_LZs_MEQ_models_nlme.R and LZs_alphaPSD_predicts_MEQ4_time_model_nlme.R first.",
      call. = FALSE
    )
  }
}

alpha_eye <- read.csv(alpha_eye_csv, stringsAsFactors = FALSE, colClasses = c(Participant = "character"))
alpha_meq <- read.csv(alpha_meq_csv, stringsAsFactors = FALSE, colClasses = c(Participant = "character"))
alpha_fixed <- read.csv(alpha_fixed_csv, stringsAsFactors = FALSE)
fixed <- read.csv(fixed_csv, stringsAsFactors = FALSE)
increment <- read.csv(increment_csv, stringsAsFactors = FALSE)

if (!requireNamespace("nlme", quietly = TRUE)) {
  stop("The nlme package is required for poster model-effect figures.", call. = FALSE)
}

epoch_levels_all <- paste0("Epoch", 0:4)
epoch_levels_dosing <- paste0("Epoch", 1:4)
epoch_labels_all <- c(
  Epoch0 = "Baseline",
  Epoch1 = "90 min",
  Epoch2 = "180 min",
  Epoch3 = "270 min",
  Epoch4 = "360 min"
)
epoch_labels_dosing <- epoch_labels_all[epoch_levels_dosing]

col_ec <- "#1f77b4"
col_eo <- "#d62728"
col_lzs <- "#d62728"
col_psd <- "#1f77b4"
col_meq <- "#2ca02c"
col_gray <- "gray35"

poster_cex_axis <- 1.45
poster_cex_lab <- 1.70
poster_cex_main <- 1.55
poster_cex_legend <- 1.30
poster_cex_note <- 1.05
poster_cex_sig <- 1.25
poster_cex_small_axis <- 1.05

theme_par <- function(mar = c(7.0, 7.8, 5.2, 12.5)) {
  par(
    bg = "white", fg = "black", col.axis = "black", col.lab = "black",
    col.main = "black", family = "sans", las = 1, bty = "l", mar = mar,
    xpd = FALSE,
    cex.lab = poster_cex_lab,
    cex.axis = poster_cex_axis,
    cex.main = poster_cex_main
  )
}

sig_label <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

format_p <- function(p) {
  vapply(p, function(one_p) {
    if (!is.finite(one_p)) return("p = NA")
    if (one_p < 0.001) return("p < .001")
    paste0("p = ", sub("^0", "", sprintf("%.3f", one_p)))
  }, character(1))
}

mean_ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) {
    return(c(N = 0, Mean = NA_real_, SE = NA_real_, CI95_low = NA_real_, CI95_high = NA_real_))
  }
  se <- if (n > 1) sd(x) / sqrt(n) else NA_real_
  ci <- if (n > 1) qt(0.975, df = n - 1) * se else NA_real_
  c(N = n, Mean = mean(x), SE = se, CI95_low = mean(x) - ci, CI95_high = mean(x) + ci)
}

summary_by <- function(data, value_col, group_cols) {
  keys <- unique(data[group_cols])
  out <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    keep <- rep(TRUE, nrow(data))
    for (col in group_cols) keep <- keep & data[[col]] == keys[[col]][i]
    stats <- mean_ci(data[[value_col]][keep])
    out[[i]] <- cbind(keys[i, , drop = FALSE], as.data.frame(as.list(stats)))
  }
  do.call(rbind, out)
}

paired_epoch_tests <- function(data, value_col) {
  rows <- list()
  for (epoch in epoch_levels_all) {
    sub <- data[data$Epoch == epoch, ]
    wide <- reshape(
      sub[, c("Participant", "Eyes", value_col)],
      idvar = "Participant",
      timevar = "Eyes",
      direction = "wide"
    )
    ec_col <- paste0(value_col, ".EC")
    eo_col <- paste0(value_col, ".EO")
    keep <- is.finite(wide[[ec_col]]) & is.finite(wide[[eo_col]])
    if (sum(keep) >= 3) {
      tt <- t.test(wide[[eo_col]][keep], wide[[ec_col]][keep], paired = TRUE)
      p <- tt$p.value
    } else {
      p <- NA_real_
    }
    rows[[length(rows) + 1]] <- data.frame(Epoch = epoch, N = sum(keep), p_value = p)
  }
  do.call(rbind, rows)
}

draw_outside_legend <- function(labels, cols, lty = 1, pch = 19, cex = poster_cex_legend) {
  old_xpd <- par("xpd")
  par(xpd = NA)
  usr <- par("usr")
  legend(
    x = usr[2] + 0.02 * diff(usr[1:2]),
    y = usr[4],
    legend = labels, col = cols,
    lty = lty, pch = pch, lwd = 2.6, bty = "n", cex = cex
  )
  par(xpd = old_xpd)
}

draw_error_bars <- function(x, low, high, col, lwd = 1.8, cap_width = 0.08) {
  segments(x, low, x, high, col = col, lwd = lwd)
  segments(x - cap_width, low, x + cap_width, low, col = col, lwd = lwd)
  segments(x - cap_width, high, x + cap_width, high, col = col, lwd = lwd)
}

format_model_p <- function(p) {
  if (!is.finite(p)) return("p = NA")
  if (p < 0.001) return("p < .001")
  paste0("p = ", sub("^0", "", sprintf("%.3f", p)))
}

fit_lme_effect <- function(data, formula_text, term, label, window_sec) {
  model_data <- data
  model_data$Participant <- factor(model_data$Participant)
  model_data$Epoch <- factor(as.character(model_data$Epoch), levels = unique(as.character(model_data$Epoch)))
  model_data$Eyes <- factor(model_data$Eyes)

  fit <- nlme::lme(
    fixed = as.formula(formula_text),
    data = model_data,
    random = ~1 | Participant,
    method = "REML",
    na.action = na.omit,
    control = nlme::lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)
  )
  tt <- summary(fit)$tTable
  if (!term %in% rownames(tt)) {
    stop("Term not found in model: ", term, call. = FALSE)
  }

  data.frame(
    WindowSec = window_sec,
    Label = label,
    Formula = formula_text,
    Term = term,
    Estimate = unname(tt[term, "Value"]),
    Std_Error = unname(tt[term, "Std.Error"]),
    DF = unname(tt[term, "DF"]),
    t_value = unname(tt[term, "t-value"]),
    p_value = unname(tt[term, "p-value"]),
    N = nrow(na.omit(model.frame(as.formula(formula_text), data = model_data))),
    N_Participants = length(unique(model_data$Participant)),
    stringsAsFactors = FALSE
  )
}

save_png <- function(file_name, width = 3000, height = 1800, res = 180) {
  png(file.path(fig_dir, file_name), width = width, height = height, res = res, bg = "white")
}

# ── 1. LZc by eye condition --------------------------------------------------

plot_lzs_eye_condition <- function(window_sec) {
  sub <- alpha_eye[alpha_eye$WindowSec == window_sec, ]
  sum_df <- summary_by(sub, "Mean_LZs", c("Epoch", "Eyes"))
  tests <- paired_epoch_tests(sub, "Mean_LZs")

  save_png(paste0("poster_lzs_eye_condition_window", window_sec, "s.png"))
  on.exit(dev.off(), add = TRUE)
  theme_par()

  ylim <- range(c(sum_df$CI95_low, sum_df$CI95_high), na.rm = TRUE)
  ylim <- ylim + c(-0.12, 0.25) * diff(ylim)
  x <- seq_along(epoch_levels_all)

  plot(
    x, rep(NA_real_, length(x)), xaxt = "n", ylim = ylim,
    xlab = "", ylab = "Mean LZc +/- 95% CI",
    main = "LZc by eye condition across baseline and dosing"
  )
  grid(col = "gray88")
  axis(1, at = x, labels = epoch_labels_all[epoch_levels_all])

  for (eye in c("EC", "EO")) {
    tmp <- sum_df[sum_df$Eyes == eye, ]
    tmp <- tmp[match(epoch_levels_all, tmp$Epoch), ]
    col <- if (eye == "EC") col_ec else col_eo
    lines(x, tmp$Mean, col = col, lwd = 2.8)
    points(x, tmp$Mean, col = col, pch = 19, cex = 1.35)
    draw_error_bars(x, tmp$CI95_low, tmp$CI95_high, col = col, lwd = 1.8)
  }

  y_sig <- ylim[2] - 0.05 * diff(ylim)
  for (i in seq_along(epoch_levels_all)) {
    p <- tests$p_value[tests$Epoch == epoch_levels_all[i]]
    n <- tests$N[tests$Epoch == epoch_levels_all[i]]
    text(x[i], y_sig, paste0(sig_label(p), "\n", "n=", n), cex = poster_cex_sig, font = 2)
  }

  draw_outside_legend(c("Eyes closed", "Eyes open"), c(col_ec, col_eo))
  mtext("*** p < .001, ** p < .01, * p < .05, ns p >= .05", side = 1, line = 4.2, cex = poster_cex_note)
}

# ── 2. Alpha PSD by eye condition -------------------------------------------

plot_alpha_eye_condition <- function(window_sec) {
  sub <- alpha_eye[alpha_eye$WindowSec == window_sec, ]
  sum_df <- summary_by(sub, "Log10Power", c("Epoch", "Eyes"))
  tests <- paired_epoch_tests(sub, "Log10Power")

  save_png(paste0("poster_alpha_psd_eye_condition_window", window_sec, "s.png"))
  on.exit(dev.off(), add = TRUE)
  theme_par()

  ylim <- range(c(sum_df$CI95_low, sum_df$CI95_high), na.rm = TRUE)
  ylim <- ylim + c(-0.12, 0.25) * diff(ylim)
  x <- seq_along(epoch_levels_all)

  plot(
    x, rep(NA_real_, length(x)), xaxt = "n", ylim = ylim,
    xlab = "", ylab = "Alpha log10 power +/- 95% CI",
    main = "Alpha PSD by eye condition across baseline and dosing"
  )
  grid(col = "gray88")
  axis(1, at = x, labels = epoch_labels_all[epoch_levels_all])

  for (eye in c("EC", "EO")) {
    tmp <- sum_df[sum_df$Eyes == eye, ]
    tmp <- tmp[match(epoch_levels_all, tmp$Epoch), ]
    col <- if (eye == "EC") col_ec else col_eo
    lines(x, tmp$Mean, col = col, lwd = 2.8)
    points(x, tmp$Mean, col = col, pch = 19, cex = 1.35)
    draw_error_bars(x, tmp$CI95_low, tmp$CI95_high, col = col, lwd = 1.8)
  }

  y_sig <- ylim[2] - 0.05 * diff(ylim)
  for (i in seq_along(epoch_levels_all)) {
    p <- tests$p_value[tests$Epoch == epoch_levels_all[i]]
    n <- tests$N[tests$Epoch == epoch_levels_all[i]]
    text(x[i], y_sig, paste0(sig_label(p), "\n", "n=", n), cex = poster_cex_sig, font = 2)
  }

  draw_outside_legend(c("Eyes closed", "Eyes open"), c(col_ec, col_eo))
  mtext("Paired EO vs EC tests at each epoch", side = 3, line = 0.2, cex = poster_cex_note)
}

# ── 3. Baseline-to-dosing trajectory used by MEQ time models -----------------

plot_meq_time_trajectory <- function(window_sec) {
  eeg_eye <- alpha_eye[alpha_eye$WindowSec == window_sec, ]
  eeg_meq <- alpha_meq[alpha_meq$WindowSec == window_sec, ]

  eeg_epoch <- aggregate(
    cbind(z_Mean_LZs, z_Log10Power) ~ Participant + Epoch,
    data = eeg_eye,
    FUN = mean,
    na.rm = TRUE
  )

  meq_epoch <- eeg_meq[, c("Participant", "Epoch", "MEQ4")]
  meq_epoch$z_MEQ4 <- as.numeric(scale(meq_epoch$MEQ4))

  sub <- merge(eeg_epoch, meq_epoch[, c("Participant", "Epoch", "z_MEQ4")],
    by = c("Participant", "Epoch"),
    all.x = TRUE
  )

  vars <- c("z_Mean_LZs", "z_Log10Power", "z_MEQ4")
  var_labels <- c(z_Mean_LZs = "LZc", z_Log10Power = "Alpha PSD", z_MEQ4 = "MEQ4")
  var_cols <- c(z_Mean_LZs = col_lzs, z_Log10Power = col_psd, z_MEQ4 = col_meq)

  rows <- list()
  for (epoch in epoch_levels_all) {
    epoch_sub <- sub[sub$Epoch == epoch, ]
    for (var in vars) {
      stats <- mean_ci(epoch_sub[[var]])
      rows[[length(rows) + 1]] <- data.frame(
        Epoch = epoch, Variable = var, N = stats["N"], Mean = stats["Mean"],
        CI95_low = stats["CI95_low"], CI95_high = stats["CI95_high"]
      )
    }
  }
  sum_df <- do.call(rbind, rows)

  save_png(paste0("poster_meq_time_model_trajectory_window", window_sec, "s.png"))
  on.exit(dev.off(), add = TRUE)
  theme_par()

  ylim <- range(c(sum_df$CI95_low, sum_df$CI95_high), na.rm = TRUE)
  ylim <- ylim + c(-0.12, 0.18) * diff(ylim)
  x <- seq_along(epoch_levels_all)

  plot(
    x, rep(NA_real_, length(x)), xaxt = "n", ylim = ylim,
    xlab = "", ylab = "Standardized mean +/- 95% CI",
    main = "Baseline and dosing trajectories for MEQ4 time models"
  )
  grid(col = "gray88")
  axis(1, at = x, labels = epoch_labels_all[epoch_levels_all])
  abline(h = 0, lty = 2, col = "gray45")

  for (var in vars) {
    tmp <- sum_df[sum_df$Variable == var, ]
    tmp <- tmp[match(epoch_levels_all, tmp$Epoch), ]
    lines(x, tmp$Mean, col = var_cols[var], lwd = 2.8)
    points(x, tmp$Mean, col = var_cols[var], pch = 19, cex = 1.35)
    draw_error_bars(x, tmp$CI95_low, tmp$CI95_high, col = var_cols[var], lwd = 1.8)
  }

  draw_outside_legend(unname(var_labels[vars]), unname(var_cols[vars]))
  mtext("MEQ4 is unavailable at baseline, so its line begins at 90 min.", side = 1, line = 4.2, cex = poster_cex_note)
}

# ── 3b. Baseline-to-dosing LZc and alpha PSD trajectory ----------------------

plot_alpha_lzc_time_trajectory <- function(window_sec) {
  eeg_eye <- alpha_eye[alpha_eye$WindowSec == window_sec, ]

  eeg_epoch <- aggregate(
    cbind(z_Mean_LZs, z_Log10Power) ~ Participant + Epoch,
    data = eeg_eye,
    FUN = mean,
    na.rm = TRUE
  )

  vars <- c("z_Mean_LZs", "z_Log10Power")
  var_labels <- c(z_Mean_LZs = "LZc", z_Log10Power = "Alpha PSD")
  var_cols <- c(z_Mean_LZs = col_lzs, z_Log10Power = col_psd)

  rows <- list()
  for (epoch in epoch_levels_all) {
    epoch_sub <- eeg_epoch[eeg_epoch$Epoch == epoch, ]
    for (var in vars) {
      stats <- mean_ci(epoch_sub[[var]])
      rows[[length(rows) + 1]] <- data.frame(
        Epoch = epoch, Variable = var, N = stats["N"], Mean = stats["Mean"],
        CI95_low = stats["CI95_low"], CI95_high = stats["CI95_high"]
      )
    }
  }
  sum_df <- do.call(rbind, rows)

  save_png(paste0("poster_alpha_psd_lzc_time_trajectory_window", window_sec, "s.png"))
  on.exit(dev.off(), add = TRUE)
  theme_par()

  ylim <- range(c(sum_df$CI95_low, sum_df$CI95_high), na.rm = TRUE)
  ylim <- ylim + c(-0.14, 0.20) * diff(ylim)
  x <- seq_along(epoch_levels_all)

  plot(
    x, rep(NA_real_, length(x)), xaxt = "n", ylim = ylim,
    xlab = "", ylab = "Standardized mean +/- 95% CI",
    main = "Baseline and dosing trajectories for LZc and alpha PSD"
  )
  grid(col = "gray88")
  axis(1, at = x, labels = epoch_labels_all[epoch_levels_all])
  abline(h = 0, lty = 2, col = "gray45")

  for (var in vars) {
    tmp <- sum_df[sum_df$Variable == var, ]
    tmp <- tmp[match(epoch_levels_all, tmp$Epoch), ]
    lines(x, tmp$Mean, col = var_cols[var], lwd = 3.0)
    points(x, tmp$Mean, col = var_cols[var], pch = 19, cex = 1.45)
    draw_error_bars(x, tmp$CI95_low, tmp$CI95_high, col = var_cols[var], lwd = 1.9)
  }

  draw_outside_legend(unname(var_labels[vars]), unname(var_cols[vars]))
  mtext("Values are averaged across eyes within participant and epoch before plotting.", side = 1, line = 4.2, cex = poster_cex_note)
}

# ── 4. MEQ4 predictor attenuation -------------------------------------------

plot_predictor_attenuation <- function(window_sec) {
  fixed_sub <- fixed[fixed$WindowSec == window_sec, ]
  keep <- fixed_sub$Term %in% c("z_Mean_LZs", "z_Log10Power") &
    fixed_sub$ModelID %in% c("M1_simple", "M2_time_controlled") &
    fixed_sub$ModelSet %in% c("LZs_only", "AlphaPSD_only", "LZs_plus_AlphaPSD")
  coefs <- fixed_sub[keep, ]
  coefs <- coefs[!(coefs$ModelSet == "LZs_plus_AlphaPSD" & coefs$ModelID == "M1_simple"), ]
  coefs$Metric <- ifelse(coefs$Term == "z_Mean_LZs", "LZc", "Alpha PSD")
  coefs$ModelLabel <- ifelse(
    coefs$ModelID == "M1_simple", "No time\ncontrol",
    ifelse(coefs$ModelSet == "LZs_plus_AlphaPSD", "+ other EEG\npredictor\n+ time", "+ time\ncontrol")
  )
  coefs$Key <- paste(coefs$Metric, coefs$ModelLabel, sep = "|")
  order_keys <- c(
    "LZc|No time\ncontrol", "LZc|+ time\ncontrol", "LZc|+ other EEG\npredictor\n+ time",
    "Alpha PSD|No time\ncontrol", "Alpha PSD|+ time\ncontrol", "Alpha PSD|+ other EEG\npredictor\n+ time"
  )
  coefs <- coefs[match(order_keys, coefs$Key, nomatch = 0), ]
  coefs$CI95_low <- coefs$Estimate - 1.96 * coefs$Std_Error
  coefs$CI95_high <- coefs$Estimate + 1.96 * coefs$Std_Error

  save_png(paste0("poster_meq_predictor_attenuation_window", window_sec, "s.png"))
  on.exit(dev.off(), add = TRUE)
  theme_par(mar = c(10.4, 7.8, 5.2, 4.2))

  x <- seq_len(nrow(coefs))
  cols <- ifelse(coefs$Metric == "LZc", col_lzs, col_psd)
  ylim <- range(c(coefs$CI95_low, coefs$CI95_high, 0), na.rm = TRUE)
  ylim <- ylim + c(-0.12, 0.22) * diff(ylim)

  plot(
    x, coefs$Estimate, type = "n", xaxt = "n", xlim = c(0.35, 6.65), ylim = ylim,
    xlab = "", ylab = "Fixed-effect estimate predicting MEQ4",
    main = "EEG predictors of MEQ4 weaken after time control"
  )
  grid(col = "gray88")
  abline(h = 0, lty = 2, col = "gray45")
  abline(v = 3.5, lty = 3, col = "gray75")
  draw_error_bars(x, coefs$CI95_low, coefs$CI95_high, col = cols, lwd = 2.2, cap_width = 0.10)
  points(x, coefs$Estimate, col = cols, pch = 19, cex = 1.45)
  axis(1, at = x, labels = coefs$ModelLabel, las = 2, cex.axis = poster_cex_small_axis)

  y_pad <- 0.04 * diff(ylim)
  text(x, coefs$CI95_high + y_pad, vapply(coefs$p_value, sig_label, character(1)), cex = poster_cex_sig, font = 2)
  text(2, ylim[2] - 0.05 * diff(ylim), "LZc predicting MEQ4", cex = poster_cex_note, font = 2, col = col_lzs)
  text(5, ylim[2] - 0.05 * diff(ylim), "Alpha PSD predicting MEQ4", cex = poster_cex_note, font = 2, col = col_psd)
  mtext("The time-controlled models add epoch; combined models include both LZc and alpha PSD.", side = 1, line = 6.4, cex = poster_cex_note)
}

# ── 5. Incremental likelihood-ratio tests -----------------------------------

plot_incremental_lrt <- function(window_sec) {
  inc <- increment[increment$WindowSec == window_sec, ]
  full_rows <- inc[!is.na(inc$Test) & inc$Test != "", ]
  if (nrow(full_rows) == 0) return(invisible(FALSE))
  labels <- ifelse(
    grepl("alpha PSD add", full_rows$ModelSet),
    "Alpha PSD added\nto LZc + Epoch",
    "LZc added\nto alpha PSD + Epoch"
  )
  lrt <- full_rows$L.Ratio
  pvals <- full_rows$p.value

  save_png(paste0("poster_meq_incremental_lrt_window", window_sec, "s.png"), width = 2600, height = 1700)
  on.exit(dev.off(), add = TRUE)
  theme_par(mar = c(9.8, 7.8, 5.2, 3.2))

  ylim <- c(0, max(lrt, na.rm = TRUE) * 1.40)
  bp <- barplot(
    lrt, names.arg = rep("", length(labels)), ylim = ylim, col = c(col_psd, col_lzs), border = NA,
    ylab = "Likelihood-ratio statistic",
    main = "Incremental MEQ4 model tests after controlling Epoch"
  )
  grid(nx = NA, ny = NULL, col = "gray88")
  old_xpd <- par("xpd")
  par(xpd = NA)
  text(bp, -0.08 * diff(ylim), labels, cex = poster_cex_axis)
  par(xpd = old_xpd)
  text(
    bp, lrt + 0.08 * max(lrt, na.rm = TRUE),
    sprintf("%s\n%s", format_p(pvals), vapply(pvals, sig_label, character(1))),
    cex = poster_cex_sig, font = 2
  )
}

# ── 6. Alpha PSD change predicted by LZc change ------------------------------

plot_alpha_lzc_change_model <- function(window_sec) {
  sub <- alpha_eye[
    alpha_eye$WindowSec == window_sec &
      alpha_eye$Epoch %in% epoch_levels_dosing &
      is.finite(alpha_eye$Delta_Log10Power) &
      is.finite(alpha_eye$Delta_Mean_LZs),
  ]
  if (nrow(sub) == 0) return(invisible(FALSE))

  model_row <- alpha_fixed[
    alpha_fixed$WindowSec == window_sec &
      alpha_fixed$Predictors == "Delta_Mean_LZs + Eyes + Epoch" &
      alpha_fixed$Term == "Delta_Mean_LZs",
  ]
  if (nrow(model_row) > 0) {
    model_text <- sprintf(
      "Mixed model controlling eye condition and time: beta = %.2f, t(%d) = %.2f, %s",
      model_row$Estimate[1], round(model_row$DF[1]), model_row$t_value[1], format_model_p(model_row$p_value[1])
    )
  } else {
    model_text <- "Mixed model controlling eye condition and time"
  }

  fit <- lm(Delta_Log10Power ~ Delta_Mean_LZs, data = sub)
  pred_x <- seq(min(sub$Delta_Mean_LZs, na.rm = TRUE), max(sub$Delta_Mean_LZs, na.rm = TRUE), length.out = 100)
  pred_y <- predict(fit, newdata = data.frame(Delta_Mean_LZs = pred_x))

  epoch_cols <- c(
    Epoch1 = "#1f77b4",
    Epoch2 = "#2ca02c",
    Epoch3 = "#d62728",
    Epoch4 = "#9467bd"
  )

  save_png(paste0("poster_alpha_psd_lzc_change_model_window", window_sec, "s.png"), width = 3200, height = 2000)
  on.exit(dev.off(), add = TRUE)
  theme_par(mar = c(7.6, 9.2, 5.4, 13.0))

  xlim <- range(sub$Delta_Mean_LZs, na.rm = TRUE)
  ylim <- range(sub$Delta_Log10Power, na.rm = TRUE)
  xlim <- xlim + c(-0.08, 0.08) * diff(xlim)
  ylim <- ylim + c(-0.14, 0.18) * diff(ylim)

  plot(
    sub$Delta_Mean_LZs, sub$Delta_Log10Power,
    type = "n", xlim = xlim, ylim = ylim,
    xlab = "Change in LZc from baseline",
    ylab = "Change in alpha log10 power from baseline",
    main = "Alpha PSD change is inversely associated with LZc change"
  )
  grid(col = "gray88")
  abline(h = 0, v = 0, lty = 2, col = "gray45")
  lines(pred_x, pred_y, col = "black", lwd = 2.4)

  for (epoch in epoch_levels_dosing) {
    tmp <- sub[sub$Epoch == epoch, ]
    points(
      tmp$Delta_Mean_LZs, tmp$Delta_Log10Power,
      pch = ifelse(tmp$Eyes == "EO", 19, 1),
      col = epoch_cols[epoch],
      bg = epoch_cols[epoch],
      cex = 1.3,
      lwd = 2
    )
  }

  old_xpd <- par("xpd")
  par(xpd = NA)
  usr <- par("usr")
  legend(
    x = usr[2] + 0.02 * diff(usr[1:2]),
    y = usr[4],
    legend = unname(epoch_labels_dosing),
    col = unname(epoch_cols[epoch_levels_dosing]),
    pch = 19, bty = "n", cex = poster_cex_note
  )
  legend(
    x = usr[2] + 0.02 * diff(usr[1:2]),
    y = usr[4] - 0.32 * diff(usr[3:4]),
    legend = c("Eyes open", "Eyes closed"),
    col = "black", pch = c(19, 1), bty = "n", cex = poster_cex_note
  )
  par(xpd = old_xpd)

  mtext(model_text, side = 3, line = 0.2, cex = poster_cex_note)
  mtext("Baseline is the zero-change reference; plotted points are dosing epochs.", side = 1, line = 4.1, cex = poster_cex_note)
}

# ── 7. Model coefficient for PSD-LZc change association ----------------------

plot_alpha_lzc_change_effect <- function(window_sec) {
  model_row <- alpha_fixed[
    alpha_fixed$WindowSec == window_sec &
      alpha_fixed$Predictors == "Delta_Mean_LZs + Eyes + Epoch" &
      alpha_fixed$Term == "Delta_Mean_LZs",
  ]
  if (nrow(model_row) == 0) return(invisible(FALSE))

  beta <- model_row$Estimate[1]
  se <- model_row$Std_Error[1]
  df <- model_row$DF[1]
  ci_width <- qt(0.975, df = df) * se
  ci_low <- beta - ci_width
  ci_high <- beta + ci_width
  p_value <- model_row$p_value[1]

  save_png(paste0("poster_alpha_psd_lzc_change_effect_window", window_sec, "s.png"), width = 2700, height = 1600)
  on.exit(dev.off(), add = TRUE)
  theme_par(mar = c(7.2, 9.0, 5.4, 3.2))

  xlim <- range(c(ci_low, ci_high, 0), na.rm = TRUE)
  xlim <- xlim + c(-0.18, 0.30) * diff(xlim)

  plot(
    beta, 1,
    type = "n",
    xlim = xlim,
    ylim = c(0.55, 1.45),
    yaxt = "n",
    xlab = "Fixed-effect estimate predicting alpha PSD change",
    ylab = "",
    main = "LZc change predicts alpha PSD change from baseline"
  )
  grid(nx = NULL, ny = NA, col = "gray88")
  abline(v = 0, lty = 2, col = "gray45")
  axis(2, at = 1, labels = "LZc change", las = 1)

  segments(ci_low, 1, ci_high, 1, col = col_lzs, lwd = 3)
  segments(ci_low, 0.94, ci_low, 1.06, col = col_lzs, lwd = 3)
  segments(ci_high, 0.94, ci_high, 1.06, col = col_lzs, lwd = 3)
  points(beta, 1, col = col_lzs, pch = 19, cex = 1.9)

  text(beta, 1.18, sig_label(p_value), cex = poster_cex_sig, font = 2)
  text(
    beta, 0.78,
    sprintf("beta = %.2f, 95%% CI [%.2f, %.2f]\nt(%d) = %.2f, %s",
      beta, ci_low, ci_high, round(df), model_row$t_value[1], format_model_p(p_value)
    ),
    cex = poster_cex_note
  )
  mtext("Model controls for eye condition, dosing time, and participant-level repeated measures.", side = 1, line = 4.1, cex = poster_cex_note)
}

# ── 7b. PSD-LZc association before and after epoch control -------------------

plot_alpha_lzc_epoch_control_effects <- function(window_sec) {
  level_data <- alpha_eye[alpha_eye$WindowSec == window_sec, ]
  level_data$Epoch <- factor(level_data$Epoch, levels = epoch_levels_all)
  level_data$Eyes <- factor(level_data$Eyes, levels = c("EC", "EO"))

  change_data <- alpha_eye[
    alpha_eye$WindowSec == window_sec &
      alpha_eye$Session == "dosing" &
      alpha_eye$Epoch != "Epoch0",
  ]
  change_data$Epoch <- factor(change_data$Epoch, levels = epoch_levels_dosing)
  change_data$Eyes <- factor(change_data$Eyes, levels = c("EC", "EO"))

  rows <- list(
    fit_lme_effect(level_data, "z_Log10Power ~ z_Mean_LZs", "z_Mean_LZs", "No\ncontrols", window_sec),
    fit_lme_effect(level_data, "z_Log10Power ~ z_Mean_LZs + Eyes", "z_Mean_LZs", "+ eyes", window_sec),
    fit_lme_effect(level_data, "z_Log10Power ~ z_Mean_LZs + Eyes + Epoch", "z_Mean_LZs", "+ eyes\n+ time", window_sec),
    fit_lme_effect(change_data, "z_Delta_Log10Power ~ z_Delta_Mean_LZs", "z_Delta_Mean_LZs", "No\ncontrols", window_sec),
    fit_lme_effect(change_data, "z_Delta_Log10Power ~ z_Delta_Mean_LZs + Eyes", "z_Delta_Mean_LZs", "+ eyes", window_sec),
    fit_lme_effect(change_data, "z_Delta_Log10Power ~ z_Delta_Mean_LZs + Eyes + Epoch", "z_Delta_Mean_LZs", "+ eyes\n+ time", window_sec)
  )
  effects <- do.call(rbind, rows)
  effects$CI95_low <- effects$Estimate - qt(0.975, df = effects$DF) * effects$Std_Error
  effects$CI95_high <- effects$Estimate + qt(0.975, df = effects$DF) * effects$Std_Error
  effects$Family <- rep(c("Level association", "Change from baseline"), each = 3)

  write.csv(
    effects,
    file.path(output_dir, paste0("poster_alpha_psd_lzc_epoch_control_effects_window", window_sec, "s.csv")),
    row.names = FALSE
  )

  save_png(paste0("poster_alpha_psd_lzc_epoch_control_effects_window", window_sec, "s.png"), width = 3300, height = 2000)
  on.exit(dev.off(), add = TRUE)
  theme_par(mar = c(10.2, 9.2, 5.4, 3.2))

  x <- c(1, 2, 3, 4.8, 5.8, 6.8)
  xlim <- c(0.35, 7.45)
  ylim <- range(c(effects$CI95_low, effects$CI95_high, 0), na.rm = TRUE)
  ylim <- ylim + c(-0.14, 0.24) * diff(ylim)

  plot(
    x, effects$Estimate,
    type = "n", xaxt = "n", xlim = xlim, ylim = ylim,
    xlab = "", ylab = "Effect of LZc predicting alpha PSD (standardized beta)",
    main = "LZc-alpha PSD association remains significant after time control"
  )
  grid(nx = NA, ny = NULL, col = "gray88")
  abline(h = 0, lty = 2, col = "gray45")
  abline(v = 3.9, lty = 3, col = "gray75")

  cols <- ifelse(effects$Family == "Level association", col_lzs, col_psd)
  for (i in seq_len(nrow(effects))) {
    segments(x[i], effects$CI95_low[i], x[i], effects$CI95_high[i], col = cols[i], lwd = 3)
    segments(x[i] - 0.08, effects$CI95_low[i], x[i] + 0.08, effects$CI95_low[i], col = cols[i], lwd = 3)
    segments(x[i] - 0.08, effects$CI95_high[i], x[i] + 0.08, effects$CI95_high[i], col = cols[i], lwd = 3)
    points(x[i], effects$Estimate[i], col = cols[i], pch = 19, cex = 1.65)
    text(x[i], effects$CI95_high[i] + 0.055 * diff(ylim), sig_label(effects$p_value[i]), cex = poster_cex_sig, font = 2)
  }

  axis(1, at = x, labels = effects$Label, las = 2, cex.axis = poster_cex_small_axis)
  text(mean(x[1:3]), ylim[2] - 0.05 * diff(ylim), "All baseline + dosing observations", cex = poster_cex_note, font = 2, col = col_lzs)
  text(mean(x[4:6]), ylim[2] - 0.05 * diff(ylim), "Baseline-to-dosing change scores", cex = poster_cex_note, font = 2, col = col_psd)
  mtext("Both sections model alpha PSD as the outcome and LZc as the predictor.", side = 1, line = 4.6, cex = poster_cex_note)
  mtext("*** p < .001, ** p < .01, * p < .05, ns p >= .05", side = 1, line = 5.7, cex = poster_cex_note)
}

# ── 8. Single poster summary of final mixed-effects conclusions --------------

plot_mixed_effects_conclusion_summary <- function(window_sec) {
  alpha_row <- alpha_fixed[
    alpha_fixed$WindowSec == window_sec &
      alpha_fixed$Predictors == "Delta_Mean_LZs + Eyes + Epoch" &
      alpha_fixed$Term == "Delta_Mean_LZs",
  ]
  lzc_meq_row <- fixed[
    fixed$WindowSec == window_sec &
      fixed$ModelSet == "LZs_only" &
      fixed$ModelID == "M2_time_controlled" &
      fixed$Term == "z_Mean_LZs",
  ]
  psd_meq_row <- fixed[
    fixed$WindowSec == window_sec &
      fixed$ModelSet == "AlphaPSD_only" &
      fixed$ModelID == "M2_time_controlled" &
      fixed$Term == "z_Log10Power",
  ]
  if (nrow(alpha_row) == 0 || nrow(lzc_meq_row) == 0 || nrow(psd_meq_row) == 0) {
    return(invisible(FALSE))
  }

  rows <- data.frame(
    Label = c(
      "LZc change -> alpha PSD change",
      "LZc -> MEQ4, controlling time",
      "Alpha PSD -> MEQ4, controlling time"
    ),
    Estimate = c(alpha_row$Estimate[1], lzc_meq_row$Estimate[1], psd_meq_row$Estimate[1]),
    DF = c(alpha_row$DF[1], lzc_meq_row$DF[1], psd_meq_row$DF[1]),
    T = c(alpha_row$t_value[1], lzc_meq_row$t_value[1], psd_meq_row$t_value[1]),
    P = c(alpha_row$p_value[1], lzc_meq_row$p_value[1], psd_meq_row$p_value[1]),
    stringsAsFactors = FALSE
  )
  rows$Evidence <- -log10(pmax(rows$P, .Machine$double.xmin))
  rows$Direction <- ifelse(rows$Estimate < 0, "negative", "positive")
  rows$Significant <- rows$P < 0.05

  save_png(paste0("poster_mixed_effects_conclusion_summary_window", window_sec, "s.png"), width = 3800, height = 2100)
  on.exit(dev.off(), add = TRUE)
  theme_par(mar = c(7.4, 15.8, 5.6, 7.0))

  y <- seq_len(nrow(rows))
  xlim <- c(0, max(rows$Evidence, -log10(0.05), na.rm = TRUE) * 2.15)
  bar_cols <- ifelse(rows$Significant, col_lzs, "gray70")

  plot(
    rows$Evidence, y,
    type = "n",
    xlim = xlim,
    ylim = c(0.45, nrow(rows) + 0.55),
    yaxt = "n",
    xlab = "Statistical evidence (-log10 p)",
    ylab = "",
    main = "Mixed-effects models support LZc-alpha PSD coupling, not MEQ4 prediction"
  )
  grid(nx = NULL, ny = NA, col = "gray88")
  abline(v = -log10(0.05), lty = 2, col = "gray45")
  axis(2, at = y, labels = rows$Label, las = 1, cex.axis = poster_cex_small_axis)

  for (i in seq_len(nrow(rows))) {
    segments(0, y[i], rows$Evidence[i], y[i], col = bar_cols[i], lwd = 14, lend = "butt")
    points(rows$Evidence[i], y[i], pch = 19, col = ifelse(rows$Significant[i], col_lzs, "gray35"), cex = 1.4)
    text(
      rows$Evidence[i] + 0.03 * diff(xlim), y[i],
      sprintf("%s  beta = %.2f, t(%d) = %.2f, %s",
        sig_label(rows$P[i]), rows$Estimate[i], round(rows$DF[i]), rows$T[i], format_model_p(rows$P[i])
      ),
      adj = 0,
      cex = poster_cex_note,
      font = ifelse(rows$Significant[i], 2, 1)
    )
  }

  text(-log10(0.05), 0.55, "p = .05", pos = 4, cex = poster_cex_note, col = "gray35")
  mtext("Only the LZc-alpha PSD model clears the p < .05 threshold; MEQ4 associations do not after time control.", side = 1, line = 4.5, cex = poster_cex_note)
}

for (window_sec in sort(unique(alpha_eye$WindowSec))) {
  plot_lzs_eye_condition(window_sec)
  plot_alpha_eye_condition(window_sec)
  plot_meq_time_trajectory(window_sec)
  plot_alpha_lzc_time_trajectory(window_sec)
  plot_predictor_attenuation(window_sec)
  plot_incremental_lrt(window_sec)
  plot_alpha_lzc_change_model(window_sec)
  plot_alpha_lzc_change_effect(window_sec)
  plot_alpha_lzc_epoch_control_effects(window_sec)
  plot_mixed_effects_conclusion_summary(window_sec)
}

legacy_alpha_file <- file.path(fig_dir, "poster_alpha_psd_eye_condition_window2s.png")
if (file.exists(legacy_alpha_file)) {
  file.copy(
    legacy_alpha_file,
    file.path(fig_dir, "poster_alpha_psd_eye_condition.png"),
    overwrite = TRUE
  )
}

cat("Wrote poster-quality figures to:\n", fig_dir, "\n\n", sep = "")
cat("Recommended poster candidates:\n")
cat("  poster_alpha_psd_eye_condition_window2s.png\n")
cat("  poster_alpha_psd_eye_condition_window10s.png\n")
cat("  poster_lzs_eye_condition_window2s.png\n")
cat("  poster_alpha_psd_lzc_time_trajectory_window2s.png\n")
cat("  poster_meq_time_model_trajectory_window2s.png\n")
cat("  poster_meq_predictor_attenuation_window2s.png\n")
cat("  poster_alpha_psd_lzc_epoch_control_effects_window2s.png\n")
cat("  poster_alpha_psd_lzc_change_model_window2s.png\n")
cat("  poster_alpha_psd_lzc_change_effect_window2s.png\n")
cat("  poster_mixed_effects_conclusion_summary_window2s.png\n")
