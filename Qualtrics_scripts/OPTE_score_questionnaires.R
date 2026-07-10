# OPTE_score_questionnaires.R

library(dplyr)

# Setup: load scoring functions
source("QA_Library.R")

# Setup: check clean files
list.files("OPTE_Qualtrics_clean")

# Setup: create scored folder
dir.create("OPTE_Qualtrics_scored", showWarnings = FALSE)


# Read cleaned files
initialprep_clean <- read.csv(
  "OPTE_Qualtrics_clean/01_OPTE_InitialPrep_CLEAN.csv",
  check.names = FALSE
)

baseline_clean <- read.csv(
  "OPTE_Qualtrics_clean/03_OPTE_Baseline_CLEAN.csv",
  check.names = FALSE
)

sessionpre_clean <- read.csv(
  "OPTE_Qualtrics_clean/04A_OPTE_SessionPRE_CLEAN.csv",
  check.names = FALSE
)

sessionpost_clean <- read.csv(
  "OPTE_Qualtrics_clean/04B_OPTE_SessionPOST_CLEAN.csv",
  check.names = FALSE
)

day1_clean <- read.csv(
  "OPTE_Qualtrics_clean/05_OPTE_1DayPostS1_CLEAN.csv",
  check.names = FALSE
)

week1_clean <- read.csv(
  "OPTE_Qualtrics_clean/06_OPTE_1WeekPostS1_CLEAN.csv",
  check.names = FALSE
)

month1_clean <- read.csv(
  "OPTE_Qualtrics_clean/07_OPTE_1MonthPostS1_CLEAN.csv",
  check.names = FALSE
)


# Check cleaned file sizes
dim(initialprep_clean)
dim(baseline_clean)
dim(sessionpre_clean)
dim(sessionpost_clean)
dim(day1_clean)
dim(week1_clean)
dim(month1_clean)


# Setup: helper function for standard questionnaire scoring
score_standard_questionnaires <- function(data) {
  
  stai_cols <- paste0("STAI_", 1:40)
  panas_cols <- paste0("PANAS_", 1:60)
  poms_cols <- paste0("POMS_", 1:65)
  bai_cols <- paste0("BAI_", 1:20)
  dass_cols <- paste0("DASS_", 1:21)
  dpes_cols <- paste0("DPES_", 1:38)
  bis11_cols <- paste0("BIS11_", 1:30)
  
  data[stai_cols] <- lapply(data[stai_cols], as.numeric)
  data[panas_cols] <- lapply(data[panas_cols], as.numeric)
  data[poms_cols] <- lapply(data[poms_cols], as.numeric)
  data[bai_cols] <- lapply(data[bai_cols], as.numeric)
  data[dass_cols] <- lapply(data[dass_cols], as.numeric)
  data[dpes_cols] <- lapply(data[dpes_cols], as.numeric)
  data[bis11_cols] <- lapply(data[bis11_cols], as.numeric)
  
  stai_output <- stai.scoring(data)
  panas_output <- panas.scoring(data)
  poms_output <- poms.scoring(data)
  bai_output <- bai.scoring(data)
  dass_output <- dass.scoring(data)
  dpes_output <- dpes.scoring(data)
  bis11_output <- bis11.scoring(data)
  
  scored <- stai_output %>%
    left_join(panas_output, by = c("ID", "date")) %>%
    left_join(poms_output, by = c("ID", "date")) %>%
    left_join(bai_output, by = c("ID", "date")) %>%
    left_join(dass_output, by = c("ID", "date")) %>%
    left_join(dpes_output, by = c("ID", "date")) %>%
    left_join(bis11_output, by = c("ID", "date"))
  
  # Keep VEQ items if this file has them
  if ("VEQ_1" %in% names(data)) {
    veq_data <- data %>%
      select(ID, date, starts_with("VEQ_"))
    
    scored <- scored %>%
      left_join(veq_data, by = c("ID", "date"))
  }
  
  return(scored)
}


# Setup: helper function for OPTE MEQ scoring
# OPTE MEQ items are already coded 0-5, so do not subtract 1
score_opte_meq <- function(data) {
  
  meq_cols <- paste0("MEQ_", 1:30)
  data[meq_cols] <- lapply(data[meq_cols], as.numeric)
  
  raw_dat <- data[, meq_cols]
  
  vars.meqMystical <- paste0("MEQ_", c(4, 5, 6, 9, 14, 15, 16, 18, 20, 21, 23, 24, 25, 26, 28))
  vars.meqPositive <- paste0("MEQ_", c(2, 8, 12, 17, 27, 30))
  vars.meqTranscendence <- paste0("MEQ_", c(1, 7, 11, 13, 19, 22))
  vars.meqIneffibility <- paste0("MEQ_", c(3, 10, 29))
  vars.meqTotal <- paste0("MEQ_", 1:30)
  
  meq.total <- rowSums(data[, vars.meqTotal]) / length(vars.meqTotal)
  meq.mystical <- rowSums(data[, vars.meqMystical]) / length(vars.meqMystical)
  meq.positive <- rowSums(data[, vars.meqPositive]) / length(vars.meqPositive)
  meq.transendence <- rowSums(data[, vars.meqTranscendence]) / length(vars.meqTranscendence)
  meq.ineffibility <- rowSums(data[, vars.meqIneffibility]) / length(vars.meqIneffibility)
  
  meq.total.pctmax <- meq.total / 5
  meq.mystical.pctmax <- meq.mystical / 5
  meq.positive.pctmax <- meq.positive / 5
  meq.transendence.pctmax <- meq.transendence / 5
  meq.ineffibility.pctmax <- meq.ineffibility / 5
  
  meq.complete <- meq.mystical.pctmax >= 0.6 &
    meq.positive.pctmax >= 0.6 &
    meq.transendence.pctmax >= 0.6 &
    meq.ineffibility.pctmax >= 0.6
  
  output <- cbind(
    data[, c("ID", "date")],
    raw_dat,
    meq.total,
    meq.mystical,
    meq.positive,
    meq.transendence,
    meq.ineffibility,
    meq.total.pctmax,
    meq.mystical.pctmax,
    meq.positive.pctmax,
    meq.transendence.pctmax,
    meq.ineffibility.pctmax,
    meq.complete
  )
  
  output <- as.data.frame(output)
  return(output)
}


# 01 Initial Prep: score questionnaires
initialprep_scored <- score_standard_questionnaires(initialprep_clean)

# 01 Initial Prep: check scored data
dim(initialprep_scored)
names(initialprep_scored)

# 01 Initial Prep: save scored file
write.csv(
  initialprep_scored,
  "OPTE_Qualtrics_scored/01_OPTE_InitialPrep_SCORED.csv",
  row.names = FALSE
)

# 01 Initial Prep: confirm scored file saved
file.exists("OPTE_Qualtrics_scored/01_OPTE_InitialPrep_SCORED.csv")



# 03 Baseline: keep cleaned VEQ data
baseline_scored <- baseline_clean

# 03 Baseline: check data
dim(baseline_scored)
names(baseline_scored)

# 03 Baseline: save file
write.csv(
  baseline_scored,
  "OPTE_Qualtrics_scored/03_OPTE_Baseline_SCORED.csv",
  row.names = FALSE
)

# 03 Baseline: confirm file saved
file.exists("OPTE_Qualtrics_scored/03_OPTE_Baseline_SCORED.csv")



# 04A Session PRE: keep cleaned PreSessQA data
sessionpre_scored <- sessionpre_clean

# 04A Session PRE: check data
dim(sessionpre_scored)
names(sessionpre_scored)

# 04A Session PRE: save file
write.csv(
  sessionpre_scored,
  "OPTE_Qualtrics_scored/04A_OPTE_SessionPRE_SCORED.csv",
  row.names = FALSE
)

# 04A Session PRE: confirm file saved
file.exists("OPTE_Qualtrics_scored/04A_OPTE_SessionPRE_SCORED.csv")



# 04B Session POST: check MEQ response range
meq_cols <- paste0("MEQ_", 1:30)

range(
  as.numeric(unlist(sessionpost_clean[meq_cols])),
  na.rm = TRUE
)

# 04B Session POST: score MEQ
# OPTE MEQ items are already coded 0-5, so use OPTE-specific MEQ scoring
sessionpost_meq <- score_opte_meq(sessionpost_clean)

# 04B Session POST: score CEQ
ceq_cols <- paste0("CEQ_", 1:26)
sessionpost_clean[ceq_cols] <- lapply(sessionpost_clean[ceq_cols], as.numeric)
sessionpost_ceq <- ceq.scoring(sessionpost_clean)

# 04B Session POST: score PsychInsight
psychinsight_cols <- paste0("PsychInsight_", 1:23)
sessionpost_clean[psychinsight_cols] <- lapply(sessionpost_clean[psychinsight_cols], as.numeric)
sessionpost_psychinsight <- psychinsight.scoring(sessionpost_clean)

# 04B Session POST: score 5DASC
fdasc_cols <- paste0("5DASC_", 1:43)
sessionpost_clean[fdasc_cols] <- lapply(sessionpost_clean[fdasc_cols], as.numeric)
sessionpost_fdasc <- fdasc.scoring(sessionpost_clean)

# 04B Session POST: combine scored outputs
sessionpost_scored <- sessionpost_meq %>%
  left_join(sessionpost_ceq, by = c("ID", "date")) %>%
  left_join(sessionpost_psychinsight, by = c("ID", "date")) %>%
  left_join(sessionpost_fdasc, by = c("ID", "date"))

# 04B Session POST: check scored data
dim(sessionpost_scored)
names(sessionpost_scored)

# 04B Session POST: save scored file
write.csv(
  sessionpost_scored,
  "OPTE_Qualtrics_scored/04B_OPTE_SessionPOST_SCORED.csv",
  row.names = FALSE
)

# 04B Session POST: confirm scored file saved
file.exists("OPTE_Qualtrics_scored/04B_OPTE_SessionPOST_SCORED.csv")



# 05 1-Day PostS1: score questionnaires
day1_scored <- score_standard_questionnaires(day1_clean)

# 05 1-Day PostS1: check scored data
dim(day1_scored)
names(day1_scored)

# 05 1-Day PostS1: save scored file
write.csv(
  day1_scored,
  "OPTE_Qualtrics_scored/05_OPTE_1DayPostS1_SCORED.csv",
  row.names = FALSE
)

# 05 1-Day PostS1: confirm scored file saved
file.exists("OPTE_Qualtrics_scored/05_OPTE_1DayPostS1_SCORED.csv")



# 06 1-Week PostS1: score questionnaires
week1_scored <- score_standard_questionnaires(week1_clean)

# 06 1-Week PostS1: check scored data
dim(week1_scored)
names(week1_scored)

# 06 1-Week PostS1: save scored file
write.csv(
  week1_scored,
  "OPTE_Qualtrics_scored/06_OPTE_1WeekPostS1_SCORED.csv",
  row.names = FALSE
)

# 06 1-Week PostS1: confirm scored file saved
file.exists("OPTE_Qualtrics_scored/06_OPTE_1WeekPostS1_SCORED.csv")



# 07 1-Month PostS1: score questionnaires
month1_scored <- score_standard_questionnaires(month1_clean)

# 07 1-Month PostS1: check scored data
dim(month1_scored)
names(month1_scored)

# 07 1-Month PostS1: save scored file
write.csv(
  month1_scored,
  "OPTE_Qualtrics_scored/07_OPTE_1MonthPostS1_SCORED.csv",
  row.names = FALSE
)

# 07 1-Month PostS1: confirm scored file saved
file.exists("OPTE_Qualtrics_scored/07_OPTE_1MonthPostS1_SCORED.csv")



# Setup: check scored files
list.files("OPTE_Qualtrics_scored")