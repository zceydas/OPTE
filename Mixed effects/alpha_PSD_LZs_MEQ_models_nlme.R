# =============================================================================
# OPTE alpha PSD, LZs, and MEQ4 linkage models using base R + nlme.
#
# This script focuses the broad PSD/LZs dataset onto alpha and asks:
#   1. Do alpha PSD and LZs line up across baseline/dosing/eyes?
#   2. Do alpha PSD changes and LZs changes line up during dosing?
#   3. Do alpha PSD and LZs relate to concurrent MEQ4 scores?
#
# Notes:
#   - "LZs" here is the Schartner-style channel/window Lempel-Ziv measure from
#     the V7 long output, averaged into Participant x Eyes x Epoch x WindowSec.
#   - MEQ4 has no eye condition, so MEQ models use values averaged across EC/EO.
#   - MEQ time mapping follows the existing project script:
#       T4 -> Epoch1, T6 -> Epoch2, T8 -> Epoch3, T9 -> Epoch4
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
fig_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

combined_csv <- file.path(output_dir, "analysis_ready_psd_band_lzs_by_epoch_window.csv")
env_combined_csv <- Sys.getenv("OPTE_PSD_LZS_READY_CSV", unset = "")
if (nzchar(env_combined_csv)) {
  combined_csv <- env_combined_csv
}

if (!file.exists(combined_csv)) {
  stop(
    "Analysis-ready PSD + LZs CSV not found:\n",
    combined_csv,
    "\nRun build_PSD_LZs_mixed_dataset_baseR.R first.",
    call. = FALSE
  )
}

cat("\nUsing PSD + LZs CSV:\n", combined_csv, "\n", sep = "")
cat("Writing alpha PSD/LZs/MEQ outputs to:\n", output_dir, "\n\n", sep = "")

combined <- read.csv(combined_csv, stringsAsFactors = FALSE, colClasses = c(Participant = "character"))

required_cols <- c(
  "Participant", "Session", "Eyes", "Epoch", "Band", "WindowSec",
  "Log10Power", "RelativePower", "Delta_Log10Power", "Delta_RelativePower",
  "z_Log10Power", "z_RelativePower", "Mean_LZs", "Mean_LZsN",
  "z_Mean_LZs", "z_Mean_LZsN"
)
missing_cols <- setdiff(required_cols, names(combined))
if (length(missing_cols) > 0) {
  stop("Combined table is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

clean_participant <- function(x) sprintf("%03d", as.integer(as.character(x)))
epoch_order <- paste0("Epoch", 0:4)

combined$Participant <- clean_participant(combined$Participant)
combined$Session <- tolower(trimws(as.character(combined$Session)))
combined$Eyes <- toupper(trimws(as.character(combined$Eyes)))
combined$Epoch <- trimws(as.character(combined$Epoch))
combined$Band <- tolower(trimws(as.character(combined$Band)))
combined$WindowSec <- as.integer(combined$WindowSec)

alpha <- combined[combined$Band == "alpha", ]
alpha <- alpha[order(alpha$Participant, alpha$Eyes, alpha$Epoch, alpha$WindowSec), ]

baseline_lzs <- alpha[alpha$Session == "baseline" & alpha$Epoch == "Epoch0",
                      c("Participant", "Eyes", "WindowSec", "Mean_LZs", "Mean_LZsN")]
names(baseline_lzs)[names(baseline_lzs) == "Mean_LZs"] <- "Baseline_Mean_LZs"
names(baseline_lzs)[names(baseline_lzs) == "Mean_LZsN"] <- "Baseline_Mean_LZsN"

alpha <- merge(alpha, baseline_lzs, by = c("Participant", "Eyes", "WindowSec"), all.x = TRUE, sort = FALSE)
alpha$Delta_Mean_LZs <- alpha$Mean_LZs - alpha$Baseline_Mean_LZs
alpha$Delta_Mean_LZsN <- alpha$Mean_LZsN - alpha$Baseline_Mean_LZsN

z_by_group <- function(df, value_col, group_cols, out_col) {
  df[[out_col]] <- NA_real_
  key <- interaction(df[group_cols], drop = TRUE, lex.order = TRUE)
  for (level in levels(key)) {
    idx <- which(key == level)
    values <- df[[value_col]][idx]
    if (sum(is.finite(values)) >= 2 && is.finite(sd(values, na.rm = TRUE)) && sd(values, na.rm = TRUE) > 0) {
      df[[out_col]][idx] <- as.numeric(scale(values))
    }
  }
  df
}

alpha <- z_by_group(alpha, "Delta_Mean_LZs", c("Participant", "WindowSec"), "z_Delta_Mean_LZs")
alpha <- z_by_group(alpha, "Delta_Mean_LZsN", c("Participant", "WindowSec"), "z_Delta_Mean_LZsN")

write.csv(alpha, file.path(output_dir, "alpha_psd_lzs_by_eye_window.csv"), row.names = FALSE)

mean_aggregate <- function(df, value_cols, key_cols) {
  out <- NULL
  for (value_col in value_cols) {
    form <- as.formula(paste(value_col, "~", paste(key_cols, collapse = " + ")))
    tmp <- aggregate(form, data = df, FUN = function(x) mean(x, na.rm = TRUE))
    names(tmp)[names(tmp) == value_col] <- value_col
    if (is.null(out)) {
      out <- tmp
    } else {
      out <- merge(out, tmp, by = key_cols, all = TRUE, sort = FALSE)
    }
  }
  out
}

alpha_epoch_avg <- mean_aggregate(
  alpha,
  c(
    "Log10Power", "RelativePower", "Delta_Log10Power", "Delta_RelativePower",
    "z_Log10Power", "z_RelativePower", "Mean_LZs", "Mean_LZsN",
    "z_Mean_LZs", "z_Mean_LZsN", "Delta_Mean_LZs", "Delta_Mean_LZsN",
    "z_Delta_Mean_LZs", "z_Delta_Mean_LZsN"
  ),
  c("Participant", "Session", "Epoch", "WindowSec")
)
alpha_epoch_avg <- alpha_epoch_avg[order(alpha_epoch_avg$Participant, alpha_epoch_avg$Epoch, alpha_epoch_avg$WindowSec), ]
write.csv(alpha_epoch_avg, file.path(output_dir, "alpha_psd_lzs_epoch_avg_for_meq.csv"), row.names = FALSE)

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

alpha_meq <- merge(
  alpha_epoch_avg[alpha_epoch_avg$Session == "dosing" & alpha_epoch_avg$Epoch != "Epoch0", ],
  meq_long,
  by = c("Participant", "Epoch"),
  all.x = TRUE,
  sort = FALSE
)
alpha_meq <- alpha_meq[is.finite(alpha_meq$MEQ4), ]
alpha_meq <- alpha_meq[order(alpha_meq$Participant, alpha_meq$Epoch, alpha_meq$WindowSec), ]
write.csv(alpha_meq, file.path(output_dir, "alpha_psd_lzs_meq_epoch_avg.csv"), row.names = FALSE)

correlation_row <- function(df, x_col, y_col, label, window_sec = NA, eyes = NA, dosing_only = FALSE) {
  keep <- is.finite(df[[x_col]]) & is.finite(df[[y_col]])
  tmp <- df[keep, ]
  if (nrow(tmp) < 4) {
    return(data.frame(
      Analysis = label, X = x_col, Y = y_col, WindowSec = window_sec, Eyes = eyes,
      DosingOnly = dosing_only, N = nrow(tmp), r = NA_real_, p_value = NA_real_,
      CI95_low = NA_real_, CI95_high = NA_real_, stringsAsFactors = FALSE
    ))
  }
  ct <- cor.test(tmp[[x_col]], tmp[[y_col]], method = "pearson")
  ci <- suppressWarnings(ct$conf.int)
  data.frame(
    Analysis = label,
    X = x_col,
    Y = y_col,
    WindowSec = window_sec,
    Eyes = eyes,
    DosingOnly = dosing_only,
    N = nrow(tmp),
    r = unname(ct$estimate),
    p_value = ct$p.value,
    CI95_low = ci[1],
    CI95_high = ci[2],
    stringsAsFactors = FALSE
  )
}

cor_rows <- list()
add_cor <- function(row) {
  cor_rows[[length(cor_rows) + 1]] <<- row
}

for (window_sec in sort(unique(alpha$WindowSec))) {
  for (eye in sort(unique(alpha$Eyes))) {
    sub <- alpha[alpha$WindowSec == window_sec & alpha$Eyes == eye, ]
    add_cor(correlation_row(sub, "z_Log10Power", "z_Mean_LZs", "Alpha PSD-LZs level by eye", window_sec, eye, FALSE))
    add_cor(correlation_row(sub, "z_RelativePower", "z_Mean_LZs", "Alpha relative PSD-LZs level by eye", window_sec, eye, FALSE))
    dosing_sub <- sub[sub$Session == "dosing" & sub$Epoch != "Epoch0", ]
    add_cor(correlation_row(dosing_sub, "Delta_Log10Power", "Delta_Mean_LZs", "Alpha PSD-LZs change by eye", window_sec, eye, TRUE))
    add_cor(correlation_row(dosing_sub, "Delta_RelativePower", "Delta_Mean_LZs", "Alpha relative PSD-LZs change by eye", window_sec, eye, TRUE))
  }
}

for (window_sec in sort(unique(alpha_epoch_avg$WindowSec))) {
  sub <- alpha_epoch_avg[alpha_epoch_avg$WindowSec == window_sec, ]
  add_cor(correlation_row(sub, "z_Log10Power", "z_Mean_LZs", "Alpha PSD-LZs level averaged across eyes", window_sec, "EC+EO average", FALSE))
  dosing_sub <- sub[sub$Session == "dosing" & sub$Epoch != "Epoch0", ]
  add_cor(correlation_row(dosing_sub, "Delta_Log10Power", "Delta_Mean_LZs", "Alpha PSD-LZs change averaged across eyes", window_sec, "EC+EO average", TRUE))
}

alpha_correlations <- do.call(rbind, cor_rows)
write.csv(alpha_correlations, file.path(output_dir, "alpha_psd_lzs_correlations.csv"), row.names = FALSE)

fit_lme_safely <- function(data, formula_text) {
  data$Participant <- factor(data$Participant)
  if ("Eyes" %in% names(data)) data$Eyes <- factor(data$Eyes, levels = c("EC", "EO"))
  data$Epoch <- factor(data$Epoch, levels = epoch_order)
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

fixed_effect_rows <- function(model, model_name, window_sec, outcome, predictors, n, n_participants) {
  if (inherits(model, "error")) {
    return(data.frame(
      Model = model_name, WindowSec = window_sec, Outcome = outcome, Predictors = predictors,
      Term = NA_character_, Estimate = NA_real_, Std_Error = NA_real_, DF = NA_real_,
      t_value = NA_real_, p_value = NA_real_, N = n, N_Participants = n_participants,
      Error = model$message, stringsAsFactors = FALSE
    ))
  }
  coefs <- as.data.frame(summary(model)$tTable)
  coefs$Term <- rownames(coefs)
  rownames(coefs) <- NULL
  data.frame(
    Model = model_name,
    WindowSec = window_sec,
    Outcome = outcome,
    Predictors = predictors,
    Term = coefs$Term,
    Estimate = coefs$Value,
    Std_Error = coefs$Std.Error,
    DF = coefs$DF,
    t_value = coefs$`t-value`,
    p_value = coefs$`p-value`,
    N = n,
    N_Participants = n_participants,
    Error = NA_character_,
    stringsAsFactors = FALSE
  )
}

model_rows <- list()
add_model <- function(row) {
  model_rows[[length(model_rows) + 1]] <<- row
}

for (window_sec in sort(unique(alpha$WindowSec))) {
  sub <- alpha[alpha$WindowSec == window_sec, ]
  sub <- sub[is.finite(sub$z_Log10Power) & is.finite(sub$z_Mean_LZs), ]
  model <- fit_lme_safely(sub, "z_Log10Power ~ z_Mean_LZs + Eyes + Epoch")
  add_model(fixed_effect_rows(model, "Alpha PSD level predicted by LZs + Eyes + Epoch", window_sec, "z_Log10Power", "z_Mean_LZs + Eyes + Epoch", nrow(sub), length(unique(sub$Participant))))

  dosing_sub <- alpha[alpha$WindowSec == window_sec & alpha$Session == "dosing" & alpha$Epoch != "Epoch0", ]
  dosing_sub <- dosing_sub[is.finite(dosing_sub$Delta_Log10Power) & is.finite(dosing_sub$Delta_Mean_LZs), ]
  dosing_sub$Epoch <- factor(dosing_sub$Epoch, levels = paste0("Epoch", 1:4))
  model <- fit_lme_safely(dosing_sub, "Delta_Log10Power ~ Delta_Mean_LZs + Eyes + Epoch")
  add_model(fixed_effect_rows(model, "Alpha PSD change predicted by LZs change + Eyes + Epoch", window_sec, "Delta_Log10Power", "Delta_Mean_LZs + Eyes + Epoch", nrow(dosing_sub), length(unique(dosing_sub$Participant))))
}

alpha_psd_lzs_models <- do.call(rbind, model_rows)
write.csv(alpha_psd_lzs_models, file.path(output_dir, "alpha_psd_lzs_nlme_fixed_effects.csv"), row.names = FALSE)

meq_model_rows <- list()
add_meq_model <- function(row) {
  meq_model_rows[[length(meq_model_rows) + 1]] <<- row
}

for (window_sec in sort(unique(alpha_meq$WindowSec))) {
  sub <- alpha_meq[alpha_meq$WindowSec == window_sec, ]
  sub$Epoch <- factor(sub$Epoch, levels = paste0("Epoch", 1:4))

  model_data <- sub[is.finite(sub$MEQ4) & is.finite(sub$z_Log10Power) & is.finite(sub$z_Mean_LZs), ]
  model <- fit_lme_safely(model_data, "MEQ4 ~ z_Log10Power + z_Mean_LZs + Epoch")
  add_meq_model(fixed_effect_rows(model, "MEQ4 predicted by alpha PSD + LZs + Epoch", window_sec, "MEQ4", "z_Log10Power + z_Mean_LZs + Epoch", nrow(model_data), length(unique(model_data$Participant))))

  model_data <- sub[is.finite(sub$MEQ4) & is.finite(sub$Delta_Log10Power) & is.finite(sub$Delta_Mean_LZs), ]
  model <- fit_lme_safely(model_data, "MEQ4 ~ Delta_Log10Power + Delta_Mean_LZs + Epoch")
  add_meq_model(fixed_effect_rows(model, "MEQ4 predicted by alpha PSD change + LZs change + Epoch", window_sec, "MEQ4", "Delta_Log10Power + Delta_Mean_LZs + Epoch", nrow(model_data), length(unique(model_data$Participant))))
}

alpha_meq_models <- do.call(rbind, meq_model_rows)
write.csv(alpha_meq_models, file.path(output_dir, "alpha_meq_psd_lzs_nlme_fixed_effects.csv"), row.names = FALSE)

plot_scatter <- function(df, x_col, y_col, file_name, main, xlab, ylab) {
  plot_df <- df[is.finite(df[[x_col]]) & is.finite(df[[y_col]]), ]
  if (nrow(plot_df) < 4) return(invisible(FALSE))
  png(file.path(fig_dir, file_name), width = 1500, height = 1200, res = 180, bg = "white")
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(bg = "white", fg = "black", col.axis = "black", col.lab = "black", col.main = "black", mar = c(5, 5, 4, 2))
  cols <- if ("Epoch" %in% names(plot_df)) {
    epoch_cols <- c(Epoch0 = "#666666", Epoch1 = "#1f77b4", Epoch2 = "#2ca02c", Epoch3 = "#d62728", Epoch4 = "#9467bd")
    epoch_cols[as.character(plot_df$Epoch)]
  } else {
    "#1f77b4"
  }
  plot(plot_df[[x_col]], plot_df[[y_col]], pch = 19, cex = 1.1, col = cols, xlab = xlab, ylab = ylab, main = main)
  grid(col = "gray88")
  abline(lm(plot_df[[y_col]] ~ plot_df[[x_col]]), lwd = 2, col = "black")
  if ("Epoch" %in% names(plot_df)) {
    legend("topright", legend = names(epoch_cols), col = epoch_cols, pch = 19, bty = "n", cex = 1.105)
  }
  invisible(TRUE)
}

for (window_sec in sort(unique(alpha_epoch_avg$WindowSec))) {
  sub <- alpha_epoch_avg[alpha_epoch_avg$WindowSec == window_sec, ]
  plot_scatter(
    sub, "z_Mean_LZs", "z_Log10Power",
    paste0("alpha_psd_lzs_level_window", window_sec, "s.png"),
    paste0("Alpha PSD and LZs level, window ", window_sec, "s"),
    "LZs z-score",
    "Alpha log10 PSD z-score"
  )
  dosing_sub <- sub[sub$Session == "dosing" & sub$Epoch != "Epoch0", ]
  plot_scatter(
    dosing_sub, "Delta_Mean_LZs", "Delta_Log10Power",
    paste0("alpha_psd_lzs_change_window", window_sec, "s.png"),
    paste0("Alpha PSD change and LZs change, window ", window_sec, "s"),
    "LZs change from baseline",
    "Alpha log10 PSD change from baseline"
  )
  meq_sub <- alpha_meq[alpha_meq$WindowSec == window_sec, ]
  plot_scatter(
    meq_sub, "z_Log10Power", "MEQ4",
    paste0("alpha_psd_meq_window", window_sec, "s.png"),
    paste0("MEQ4 and alpha PSD, window ", window_sec, "s"),
    "Alpha log10 PSD z-score",
    "MEQ4"
  )
  plot_scatter(
    meq_sub, "z_Mean_LZs", "MEQ4",
    paste0("alpha_lzs_meq_window", window_sec, "s.png"),
    paste0("MEQ4 and LZs, window ", window_sec, "s"),
    "LZs z-score",
    "MEQ4"
  )
}

cat("Wrote alpha focused data:\n")
cat("  alpha_psd_lzs_by_eye_window.csv\n")
cat("  alpha_psd_lzs_epoch_avg_for_meq.csv\n")
cat("  alpha_psd_lzs_meq_epoch_avg.csv\n\n")
cat("Wrote alpha focused results:\n")
cat("  alpha_psd_lzs_correlations.csv\n")
cat("  alpha_psd_lzs_nlme_fixed_effects.csv\n")
cat("  alpha_meq_psd_lzs_nlme_fixed_effects.csv\n\n")
cat("Wrote alpha focused figures to:\n", fig_dir, "\n", sep = "")
