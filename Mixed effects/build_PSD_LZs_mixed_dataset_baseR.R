# =============================================================================
# Build OPTE PSD + LZs mixed-model dataset using base R only.
#
# This script intentionally avoids external packages so it can run before you
# install tidyverse/lme4/lmerTest. It creates the analysis-ready CSVs needed for
# later mixed-effects models.
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
script_dir <- gsub("~+~", " ", script_dir, fixed = TRUE)
script_dir <- normalizePath(script_dir, mustWork = FALSE)
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

lzs_csv <- file.path(
  repo_root,
  "MATLAB/OPTE/data/outputs/resting_preprocessing",
  "LZc_Schartner_ProductionV7_Results/CSV/PUBLIC_long_LZs_by_channel_window.csv"
)

psd_band_csv <- file.path(
  repo_root,
  "MATLAB/OPTE/data/outputs/resting_preprocessing",
  "PSD_Results/ALL_PSD_band_power_by_channel.csv"
)

env_lzs_csv <- Sys.getenv("OPTE_LZS_CSV", unset = "")
env_psd_band_csv <- Sys.getenv("OPTE_PSD_BAND_CSV", unset = "")
if (nzchar(env_lzs_csv)) {
  lzs_csv <- env_lzs_csv
}
if (nzchar(env_psd_band_csv)) {
  psd_band_csv <- env_psd_band_csv
}

output_dir <- file.path(script_dir, "outputs", "PSD_band_LZs_mixed_models")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found:\n", path, call. = FALSE)
  }
}

clean_participant <- function(x) {
  sprintf("%03d", as.integer(as.character(x)))
}

clean_session <- function(x) {
  tolower(trimws(as.character(x)))
}

clean_eyes <- function(x) {
  toupper(trimws(as.character(x)))
}

clean_epoch <- function(x) {
  x <- trimws(as.character(x))
  lower <- tolower(x)
  out <- x
  is_epoch <- grepl("^epoch[0-9]+$", lower)
  out[is_epoch] <- paste0("Epoch", sub("^epoch", "", lower[is_epoch]))
  out
}

z_within_groups <- function(df, value_col, group_cols, out_col) {
  df[[out_col]] <- NA_real_
  group_key <- interaction(df[group_cols], drop = TRUE, lex.order = TRUE)
  for (key in levels(group_key)) {
    idx <- which(group_key == key)
    values <- df[[value_col]][idx]
    if (sum(is.finite(values)) >= 2 && is.finite(sd(values, na.rm = TRUE)) &&
        sd(values, na.rm = TRUE) > 0) {
      df[[out_col]][idx] <- as.numeric(scale(values))
    }
  }
  df
}

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

count_unique_channel <- function(df, key_cols, channel_col, out_col) {
  key <- interaction(df[key_cols], drop = TRUE, lex.order = TRUE)
  counts <- tapply(df[[channel_col]], key, function(x) length(unique(x)))
  key_df <- unique(df[key_cols])
  key_df$key_tmp <- interaction(key_df[key_cols], drop = TRUE, lex.order = TRUE)
  key_df[[out_col]] <- as.integer(counts[as.character(key_df$key_tmp)])
  key_df$key_tmp <- NULL
  key_df
}

cat("\nUsing PSD band CSV:\n", psd_band_csv, "\n", sep = "")
cat("Using LZs CSV:\n", lzs_csv, "\n", sep = "")
cat("Writing outputs to:\n", output_dir, "\n\n", sep = "")

stop_if_missing(psd_band_csv, "PSD band CSV")
stop_if_missing(lzs_csv, "LZs CSV")

# ---- PSD band table ----------------------------------------------------------

psd_raw <- read.csv(psd_band_csv, stringsAsFactors = FALSE,
                    colClasses = c(Participant = "character"))

required_psd_cols <- c(
  "Participant", "Session", "Eyes", "Epoch", "Channel", "Band",
  "AbsolutePower", "Log10Power", "RelativePower"
)
missing_psd_cols <- setdiff(required_psd_cols, names(psd_raw))
if (length(missing_psd_cols) > 0) {
  stop("PSD CSV is missing columns: ", paste(missing_psd_cols, collapse = ", "))
}

psd_raw$Participant <- clean_participant(psd_raw$Participant)
psd_raw$Session <- clean_session(psd_raw$Session)
psd_raw$Eyes <- clean_eyes(psd_raw$Eyes)
psd_raw$Epoch <- clean_epoch(psd_raw$Epoch)
psd_raw$Band <- tolower(trimws(as.character(psd_raw$Band)))

psd_keys <- c("Participant", "Session", "Eyes", "Epoch", "Band")
psd_epoch <- mean_aggregate(
  psd_raw,
  c("AbsolutePower", "Log10Power", "RelativePower"),
  psd_keys
)
psd_counts <- count_unique_channel(psd_raw, psd_keys, "Channel", "N_Channels_PSD")
psd_epoch <- merge(psd_epoch, psd_counts, by = psd_keys, all.x = TRUE, sort = FALSE)

psd_epoch <- psd_epoch[order(
  psd_epoch$Participant, psd_epoch$Session, psd_epoch$Eyes,
  psd_epoch$Epoch, psd_epoch$Band
), ]

for (col in c("AbsolutePower", "Log10Power", "RelativePower")) {
  psd_epoch <- z_within_groups(
    psd_epoch, col, c("Participant", "Band"), paste0("z_", col)
  )
}

baseline <- psd_epoch[psd_epoch$Session == "baseline" & psd_epoch$Epoch == "Epoch0",
                      c("Participant", "Eyes", "Band",
                        "AbsolutePower", "Log10Power", "RelativePower")]
names(baseline)[names(baseline) == "AbsolutePower"] <- "Baseline_AbsolutePower"
names(baseline)[names(baseline) == "Log10Power"] <- "Baseline_Log10Power"
names(baseline)[names(baseline) == "RelativePower"] <- "Baseline_RelativePower"

psd_epoch <- merge(
  psd_epoch, baseline,
  by = c("Participant", "Eyes", "Band"),
  all.x = TRUE,
  sort = FALSE
)

psd_epoch$Delta_AbsolutePower <- psd_epoch$AbsolutePower - psd_epoch$Baseline_AbsolutePower
psd_epoch$Delta_Log10Power <- psd_epoch$Log10Power - psd_epoch$Baseline_Log10Power
psd_epoch$Delta_RelativePower <- psd_epoch$RelativePower - psd_epoch$Baseline_RelativePower

for (col in c("Delta_AbsolutePower", "Delta_Log10Power", "Delta_RelativePower")) {
  psd_epoch <- z_within_groups(
    psd_epoch, col, c("Participant", "Band"), paste0("z_", col)
  )
}

write.csv(
  psd_epoch,
  file.path(output_dir, "analysis_ready_psd_band_by_epoch.csv"),
  row.names = FALSE
)

# ---- LZs table ---------------------------------------------------------------

lzs_raw <- read.csv(lzs_csv, stringsAsFactors = FALSE,
                    colClasses = c(Participant = "character"))

required_lzs_cols <- c(
  "Participant", "Session", "Eyes", "Epoch", "WindowSec", "Window",
  "Channel", "LZs", "LZsN"
)
missing_lzs_cols <- setdiff(required_lzs_cols, names(lzs_raw))
if (length(missing_lzs_cols) > 0) {
  stop("LZs CSV is missing columns: ", paste(missing_lzs_cols, collapse = ", "))
}

lzs_raw$Participant <- clean_participant(lzs_raw$Participant)
lzs_raw$Session <- clean_session(lzs_raw$Session)
lzs_raw$Eyes <- clean_eyes(lzs_raw$Eyes)
lzs_raw$Epoch <- clean_epoch(lzs_raw$Epoch)

lzs_keys <- c("Participant", "Session", "Eyes", "Epoch", "WindowSec")
lzs_epoch <- mean_aggregate(lzs_raw, c("LZs", "LZsN"), lzs_keys)
names(lzs_epoch)[names(lzs_epoch) == "LZs"] <- "Mean_LZs"
names(lzs_epoch)[names(lzs_epoch) == "LZsN"] <- "Mean_LZsN"

lzs_channel_counts <- count_unique_channel(lzs_raw, lzs_keys, "Channel", "N_Channels_LZs")
lzs_window_counts <- count_unique_channel(lzs_raw, lzs_keys, "Window", "N_Windows_LZs")
lzs_epoch <- merge(lzs_epoch, lzs_channel_counts, by = lzs_keys, all.x = TRUE, sort = FALSE)
lzs_epoch <- merge(lzs_epoch, lzs_window_counts, by = lzs_keys, all.x = TRUE, sort = FALSE)

lzs_epoch <- lzs_epoch[order(
  lzs_epoch$Participant, lzs_epoch$Session, lzs_epoch$Eyes,
  lzs_epoch$Epoch, lzs_epoch$WindowSec
), ]

lzs_epoch <- z_within_groups(lzs_epoch, "Mean_LZs", c("Participant", "WindowSec"), "z_Mean_LZs")
lzs_epoch <- z_within_groups(lzs_epoch, "Mean_LZsN", c("Participant", "WindowSec"), "z_Mean_LZsN")

write.csv(
  lzs_epoch,
  file.path(output_dir, "analysis_ready_lzs_by_epoch_window.csv"),
  row.names = FALSE
)

# ---- Combined PSD + LZs table ------------------------------------------------

combined <- merge(
  psd_epoch,
  lzs_epoch,
  by = c("Participant", "Session", "Eyes", "Epoch"),
  all.x = TRUE,
  sort = FALSE
)

combined <- combined[order(
  combined$Participant, combined$Session, combined$Eyes,
  combined$Epoch, combined$Band, combined$WindowSec
), ]

write.csv(
  combined,
  file.path(output_dir, "analysis_ready_psd_band_lzs_by_epoch_window.csv"),
  row.names = FALSE
)

cat("Wrote:\n")
cat("  analysis_ready_psd_band_by_epoch.csv\n")
cat("  analysis_ready_lzs_by_epoch_window.csv\n")
cat("  analysis_ready_psd_band_lzs_by_epoch_window.csv\n\n")

cat("Rows:\n")
cat("  PSD epoch table: ", nrow(psd_epoch), "\n", sep = "")
cat("  LZs epoch table: ", nrow(lzs_epoch), "\n", sep = "")
cat("  Combined table:  ", nrow(combined), "\n", sep = "")
cat("Participants in PSD table: ", length(unique(psd_epoch$Participant)), "\n", sep = "")
cat("Bands: ", paste(sort(unique(psd_epoch$Band)), collapse = ", "), "\n", sep = "")
