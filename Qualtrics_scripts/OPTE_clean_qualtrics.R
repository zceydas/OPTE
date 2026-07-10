# OPTE_clean_qualtrics.R

library(readr)
library(dplyr)

# Setup: check raw files
list.files("OPTE_Qualtrics_raw")

# Setup: create clean folder
dir.create("OPTE_Qualtrics_clean", showWarnings = FALSE)

# Setup: define valid participant IDs
valid_ids <- c("5", "7", "8", "11", "12", "13", "16", "17", "18", "19")

# Setup: standard questionnaire names
standard_names_276 <- c(
  "ID", "date",
  paste0("STAI_", 1:40),
  paste0("PANAS_", 1:60),
  paste0("POMS_", 1:65),
  paste0("BAI_", 1:20),
  paste0("DASS_", 1:21),
  paste0("DPES_", 1:38),
  paste0("BIS11_", 1:30)
)

standard_names_292 <- c(
  standard_names_276,
  paste0("VEQ_", 1:16)
)

sessionpost_names <- c(
  "ID", "date",
  paste0("MEQ_", 1:30),
  paste0("CEQ_", 1:26),
  paste0("PsychInsight_", 1:23),
  paste0("5DASC_", 1:43)
)


# 01 Initial Prep: read raw data
initialprep_raw <- read_csv(
  "OPTE_Qualtrics_raw/01_OPTE_InitialPrep_raw.csv",
  col_types = cols(.default = "c")
)

# 01 Initial Prep: inspect raw data
dim(initialprep_raw)
names(initialprep_raw)[1:20]

# 01 Initial Prep: inspect participant IDs and completion status
initialprep_raw %>%
  select(Status, Finished, Progress, ID, date) %>%
  print(n = Inf)

# 01 Initial Prep: clean rows
initialprep_clean <- initialprep_raw %>%
  filter(ID %in% valid_ids)

# 01 Initial Prep: clean columns
initialprep_clean <- initialprep_clean %>%
  select(
    ID, date,
    starts_with("STAI_"),
    starts_with("PANAS_"),
    starts_with("POMS_"),
    starts_with("BAI_"),
    starts_with("DASS_"),
    starts_with("DPES_"),
    starts_with("BIS11_")
  )

# 01 Initial Prep: standardize column names
names(initialprep_clean) <- standard_names_276

# 01 Initial Prep: check cleaned data
dim(initialprep_clean)
names(initialprep_clean)[1:20]

# 01 Initial Prep: save clean file
write_csv(
  initialprep_clean,
  "OPTE_Qualtrics_clean/01_OPTE_InitialPrep_CLEAN.csv"
)

# 01 Initial Prep: confirm file saved
file.exists("OPTE_Qualtrics_clean/01_OPTE_InitialPrep_CLEAN.csv")



# 03 Baseline: read raw data
baseline_raw <- read_csv(
  "OPTE_Qualtrics_raw/03_OPTE_Baseline_raw.csv",
  col_types = cols(.default = "c")
)

# 03 Baseline: inspect raw data
dim(baseline_raw)
names(baseline_raw)[1:30]

# 03 Baseline: inspect participant IDs and completion status
baseline_raw %>%
  select(Status, Finished, Progress, ID, date) %>%
  print(n = Inf)

# 03 Baseline: clean rows
baseline_clean <- baseline_raw %>%
  filter(ID %in% valid_ids)

# 03 Baseline: check IDs
baseline_clean %>%
  count(ID)

# 03 Baseline: clean columns
baseline_clean <- baseline_clean %>%
  select(
    ID, date,
    starts_with("VEQ_")
  )

# 03 Baseline: check cleaned data
dim(baseline_clean)
names(baseline_clean)

# 03 Baseline: save clean file
write_csv(
  baseline_clean,
  "OPTE_Qualtrics_clean/03_OPTE_Baseline_CLEAN.csv"
)

# 03 Baseline: confirm file saved
file.exists("OPTE_Qualtrics_clean/03_OPTE_Baseline_CLEAN.csv")



# 04A Session PRE: read raw data
sessionpre_raw <- read_csv(
  "OPTE_Qualtrics_raw/04A_OPTE_SessionPRE_raw.csv",
  col_types = cols(.default = "c")
)

# 04A Session PRE: inspect raw data
dim(sessionpre_raw)
names(sessionpre_raw)[1:30]

# 04A Session PRE: inspect participant IDs and completion status
sessionpre_raw %>%
  select(Status, Finished, Progress, ID, date) %>%
  print(n = Inf)

# 04A Session PRE: clean rows
sessionpre_clean <- sessionpre_raw %>%
  filter(ID %in% valid_ids)

# 04A Session PRE: check IDs
sessionpre_clean %>%
  count(ID)

# 04A Session PRE: clean columns
sessionpre_clean <- sessionpre_clean %>%
  select(
    ID, date,
    starts_with("PreSessQA_")
  )

# 04A Session PRE: check cleaned data
dim(sessionpre_clean)
names(sessionpre_clean)

# 04A Session PRE: save clean file
write_csv(
  sessionpre_clean,
  "OPTE_Qualtrics_clean/04A_OPTE_SessionPRE_CLEAN.csv"
)

# 04A Session PRE: confirm file saved
file.exists("OPTE_Qualtrics_clean/04A_OPTE_SessionPRE_CLEAN.csv")



# 04B Session POST: read raw data
sessionpost_raw <- read_csv(
  "OPTE_Qualtrics_raw/04B_OPTE_SessionPOST_raw.csv",
  col_types = cols(.default = "c")
)

# 04B Session POST: inspect raw data
dim(sessionpost_raw)
names(sessionpost_raw)[1:30]

# 04B Session POST: inspect participant IDs and completion status
sessionpost_raw %>%
  select(Status, Finished, Progress, ID, date) %>%
  print(n = Inf)

# 04B Session POST: clean rows
sessionpost_clean <- sessionpost_raw %>%
  filter(ID %in% valid_ids)

# 04B Session POST: check IDs
sessionpost_clean %>%
  count(ID)

# 04B Session POST: clean columns
sessionpost_clean <- sessionpost_clean %>%
  select(ID:last_col())

# 04B Session POST: standardize column names
names(sessionpost_clean) <- sessionpost_names

# 04B Session POST: check cleaned data
dim(sessionpost_clean)
names(sessionpost_clean)[1:30]
names(sessionpost_clean)[90:124]

# 04B Session POST: save clean file
write_csv(
  sessionpost_clean,
  "OPTE_Qualtrics_clean/04B_OPTE_SessionPOST_CLEAN.csv"
)

# 04B Session POST: confirm file saved
file.exists("OPTE_Qualtrics_clean/04B_OPTE_SessionPOST_CLEAN.csv")



# 05 1-Day PostS1: read raw data
day1_raw <- read_csv(
  "OPTE_Qualtrics_raw/05_OPTE_1DayPostS1_raw.csv",
  col_types = cols(.default = "c")
)

# 05 1-Day PostS1: inspect raw data
dim(day1_raw)
names(day1_raw)[1:30]

# 05 1-Day PostS1: inspect participant IDs and completion status
day1_raw %>%
  select(Status, Finished, Progress, ID, date) %>%
  print(n = Inf)

# 05 1-Day PostS1: clean rows
day1_clean <- day1_raw %>%
  filter(ID %in% valid_ids)

# 05 1-Day PostS1: check IDs
day1_clean %>%
  count(ID)

# 05 1-Day PostS1: clean columns
day1_clean <- day1_clean %>%
  select(
    ID, date,
    starts_with("STAI_"),
    starts_with("PANAS_"),
    starts_with("POMS_"),
    starts_with("BAI_"),
    starts_with("DASS_"),
    starts_with("DPES_"),
    starts_with("BIS11_")
  )

# 05 1-Day PostS1: standardize column names
names(day1_clean) <- standard_names_276

# 05 1-Day PostS1: check cleaned data
dim(day1_clean)
names(day1_clean)[1:20]

# 05 1-Day PostS1: save clean file
write_csv(
  day1_clean,
  "OPTE_Qualtrics_clean/05_OPTE_1DayPostS1_CLEAN.csv"
)

# 05 1-Day PostS1: confirm file saved
file.exists("OPTE_Qualtrics_clean/05_OPTE_1DayPostS1_CLEAN.csv")



# 06 1-Week PostS1: read raw data
week1_raw <- read_csv(
  "OPTE_Qualtrics_raw/06_OPTE_1WeekPostS1_raw.csv",
  col_types = cols(.default = "c")
)

# 06 1-Week PostS1: inspect raw data
dim(week1_raw)
names(week1_raw)[1:30]

# 06 1-Week PostS1: inspect participant IDs and completion status
week1_raw %>%
  select(Status, Finished, Progress, ID, date) %>%
  print(n = Inf)

# 06 1-Week PostS1: clean rows
week1_clean <- week1_raw %>%
  filter(ID %in% valid_ids)

# 06 1-Week PostS1: check IDs
week1_clean %>%
  count(ID)

# 06 1-Week PostS1: clean columns
week1_clean <- week1_clean %>%
  select(
    ID, date,
    starts_with("STAI_"),
    starts_with("PANAS_"),
    starts_with("POMS_"),
    starts_with("BAI_"),
    starts_with("DASS_"),
    starts_with("DPES_"),
    starts_with("BIS11_"),
    starts_with("VEQ_")
  )

# 06 1-Week PostS1: standardize column names
names(week1_clean) <- standard_names_292

# 06 1-Week PostS1: check cleaned data
dim(week1_clean)
names(week1_clean)[1:20]
names(week1_clean)[270:292]

# 06 1-Week PostS1: save clean file
write_csv(
  week1_clean,
  "OPTE_Qualtrics_clean/06_OPTE_1WeekPostS1_CLEAN.csv"
)

# 06 1-Week PostS1: confirm file saved
file.exists("OPTE_Qualtrics_clean/06_OPTE_1WeekPostS1_CLEAN.csv")



# 07 1-Month PostS1: read raw data
month1_raw <- read_csv(
  "OPTE_Qualtrics_raw/07_OPTE_1MonthPostS1_raw.csv",
  col_types = cols(.default = "c")
)

# 07 1-Month PostS1: inspect raw data
dim(month1_raw)
names(month1_raw)[1:30]

# 07 1-Month PostS1: inspect participant IDs and completion status
month1_raw %>%
  select(Status, Finished, Progress, ID, date) %>%
  print(n = Inf)

# 07 1-Month PostS1: clean rows
month1_clean <- month1_raw %>%
  filter(ID %in% valid_ids)

# 07 1-Month PostS1: check IDs before duplicate removal
month1_clean %>%
  count(ID)

# 07 1-Month PostS1: remove duplicate participant row
month1_clean <- month1_clean %>%
  arrange(ID, desc(Finished), desc(Progress)) %>%
  distinct(ID, .keep_all = TRUE)

# 07 1-Month PostS1: check IDs after duplicate removal
month1_clean %>%
  count(ID)

# 07 1-Month PostS1: clean columns
month1_clean <- month1_clean %>%
  select(
    ID, date,
    starts_with("STAI_"),
    starts_with("PANAS_"),
    starts_with("POMS_"),
    starts_with("BAI_"),
    starts_with("DASS_"),
    starts_with("DPES_"),
    starts_with("BIS11_"),
    starts_with("VEQ_")
  )

# 07 1-Month PostS1: standardize column names
names(month1_clean) <- standard_names_292

# 07 1-Month PostS1: check cleaned data
dim(month1_clean)
names(month1_clean)[1:20]
names(month1_clean)[270:292]

# 07 1-Month PostS1: save clean file
write_csv(
  month1_clean,
  "OPTE_Qualtrics_clean/07_OPTE_1MonthPostS1_CLEAN.csv"
)

# 07 1-Month PostS1: confirm file saved
file.exists("OPTE_Qualtrics_clean/07_OPTE_1MonthPostS1_CLEAN.csv")