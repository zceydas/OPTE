# OPTE_merge_scored_data.R

library(dplyr)
library(tidyr)

# Setup: check scored files
list.files("OPTE_Qualtrics_scored")

# Setup: create merged scored data folder
dir.create("OPTE_merged_scored_data", showWarnings = FALSE)


# Read scored files
initialprep_scored <- read.csv(
  "OPTE_Qualtrics_scored/01_OPTE_InitialPrep_SCORED.csv",
  check.names = FALSE
)

baseline_scored <- read.csv(
  "OPTE_Qualtrics_scored/03_OPTE_Baseline_SCORED.csv",
  check.names = FALSE
)

sessionpre_scored <- read.csv(
  "OPTE_Qualtrics_scored/04A_OPTE_SessionPRE_SCORED.csv",
  check.names = FALSE
)

sessionpost_scored <- read.csv(
  "OPTE_Qualtrics_scored/04B_OPTE_SessionPOST_SCORED.csv",
  check.names = FALSE
)

day1_scored <- read.csv(
  "OPTE_Qualtrics_scored/05_OPTE_1DayPostS1_SCORED.csv",
  check.names = FALSE
)

week1_scored <- read.csv(
  "OPTE_Qualtrics_scored/06_OPTE_1WeekPostS1_SCORED.csv",
  check.names = FALSE
)

month1_scored <- read.csv(
  "OPTE_Qualtrics_scored/07_OPTE_1MonthPostS1_SCORED.csv",
  check.names = FALSE
)


# Check scored file sizes
dim(initialprep_scored)
dim(baseline_scored)
dim(sessionpre_scored)
dim(sessionpost_scored)
dim(day1_scored)
dim(week1_scored)
dim(month1_scored)


# Add timepoint labels for LONG dataset
initialprep_long <- initialprep_scored %>%
  mutate(timepoint = "initialprep")

baseline_long <- baseline_scored %>%
  mutate(timepoint = "baseline")

sessionpre_long <- sessionpre_scored %>%
  mutate(timepoint = "sessionpre")

sessionpost_long <- sessionpost_scored %>%
  mutate(timepoint = "sessionpost")

day1_long <- day1_scored %>%
  mutate(timepoint = "1daypostS1")

week1_long <- week1_scored %>%
  mutate(timepoint = "1weekpostS1")

month1_long <- month1_scored %>%
  mutate(timepoint = "1monthpostS1")


# Create LONG dataset
opte_long <- bind_rows(
  initialprep_long,
  baseline_long,
  sessionpre_long,
  sessionpost_long,
  day1_long,
  week1_long,
  month1_long
) %>%
  select(ID, timepoint, everything())

# Check LONG dataset
dim(opte_long)
table(opte_long$timepoint)

# Save LONG dataset
write.csv(
  opte_long,
  "OPTE_merged_scored_data/OPTE_AllScoredData_LONG.csv",
  row.names = FALSE
)

# Confirm LONG file saved
file.exists("OPTE_merged_scored_data/OPTE_AllScoredData_LONG.csv")


# Helper function for WIDE dataset
add_prefix_to_columns <- function(data, prefix) {
  data %>%
    rename_with(
      .fn = function(x) paste0(prefix, "_", x),
      .cols = -ID
    )
}


# Create WIDE pieces
initialprep_wide <- add_prefix_to_columns(initialprep_scored, "initialprep")
baseline_wide <- add_prefix_to_columns(baseline_scored, "baseline")
sessionpre_wide <- add_prefix_to_columns(sessionpre_scored, "sessionpre")
sessionpost_wide <- add_prefix_to_columns(sessionpost_scored, "sessionpost")
day1_wide <- add_prefix_to_columns(day1_scored, "1daypostS1")
week1_wide <- add_prefix_to_columns(week1_scored, "1weekpostS1")
month1_wide <- add_prefix_to_columns(month1_scored, "1monthpostS1")


# Create WIDE dataset
opte_wide <- initialprep_wide %>%
  full_join(baseline_wide, by = "ID") %>%
  full_join(sessionpre_wide, by = "ID") %>%
  full_join(sessionpost_wide, by = "ID") %>%
  full_join(day1_wide, by = "ID") %>%
  full_join(week1_wide, by = "ID") %>%
  full_join(month1_wide, by = "ID") %>%
  arrange(as.numeric(ID))

# Check WIDE dataset
dim(opte_wide)
opte_wide$ID

# Save WIDE dataset
write.csv(
  opte_wide,
  "OPTE_merged_scored_data/OPTE_AllScoredData_WIDE.csv",
  row.names = FALSE
)

# Confirm WIDE file saved
file.exists("OPTE_merged_scored_data/OPTE_AllScoredData_WIDE.csv")


# Check merged scored data files
list.files("OPTE_merged_scored_data")