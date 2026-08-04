## Go/No-Go: descriptive stats, mixed-effects models, pairwise comparisons vs.
## Baseline, and Bayes Factors, for each outcome measure across time point.
##
## Random effect: participant (random intercept)
## Model form: outcome ~ period + (1 | participant)
##
## Data note: FA / zFA already use a log-linear correction for sessions with
## zero false alarms (no NaNs present in this file) -- those "no false alarm"
## sessions are kept in the analysis exactly as computed, per instructions.

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)
library(BayesFactor)

# ── Load & prepare ───────────────────────────────────────────────────────────
df <- read.csv("/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults/gonogo_metrics_by_period.csv")

period_levels <- c("Baseline TMSEEG", "Dosing Session", "1-week follow-up",
                    "2-week follow-up", "1-month follow-up")

df <- df %>%
  mutate(
    participant = factor(participant),
    period      = factor(period, levels = period_levels)
  )

stopifnot(sum(is.na(df)) == 0)

outcomes <- c("H", "FA", "zH", "zFA", "dprime", "criterion", "mean_rt")

# ── Helpers ──────────────────────────────────────────────────────────────────
banner <- paste0(rep("=", 70), collapse = "")

add_cohens_d <- function(contrast_obj, model) {
  tbl <- as.data.frame(contrast_obj)
  tbl$d <- tbl$estimate / sigma(model)
  tbl
}

paired_bf <- function(wide, s1, s2, label) {
  d <- wide[[s2]] - wide[[s1]]
  d <- d[!is.na(d)]
  if (length(d) < 3 || sd(d) == 0) {
    return(data.frame(Contrast = label, n_pairs = length(d), mean_diff = NA,
                       dz = NA, BF10 = NA, BF01 = NA))
  }
  dz   <- mean(d) / sd(d)
  bf   <- ttestBF(x = d, mu = 0)
  bf10 <- exp(bf@bayesFactor$bf)
  data.frame(Contrast = label, n_pairs = length(d), mean_diff = round(mean(d), 4),
             dz = round(dz, 3), BF10 = round(bf10, 3), BF01 = round(1 / bf10, 3))
}

# ── 1. Descriptive statistics ────────────────────────────────────────────────
desc_list <- list()
for (v in outcomes) {
  desc_list[[v]] <- df %>%
    group_by(period) %>%
    summarise(
      outcome = v,
      n    = sum(!is.na(.data[[v]])),
      mean = mean(.data[[v]], na.rm = TRUE),
      sd   = sd(.data[[v]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    relocate(outcome)
}
descriptives <- bind_rows(desc_list) %>%
  mutate(across(c(mean, sd), ~ round(.x, 4)))

write.csv(descriptives,
          "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults/gonogo_descriptives.csv",
          row.names = FALSE)

cat("\n", banner, "\nDESCRIPTIVE STATISTICS (mean / sd / n by period)\n", banner, "\n", sep = "")
print(descriptives, n = Inf)

# ── 2. Mixed models, pairwise vs Baseline, all-pairwise, Bayes Factors ─────
omnibus_list <- list()
pairwise_list <- list()
all_pairwise_list <- list()
bf_list <- list()

for (v in outcomes) {

  cat("\n\n", banner, "\nOUTCOME: ", v, "\n", banner, "\n", sep = "")

  form <- as.formula(paste0(v, " ~ period + (1 | participant)"))
  model <- lmer(form, data = df, REML = TRUE, na.action = na.omit)

  r2 <- r.squaredGLMM(model)
  cat(sprintf("\nR2m = %.3f | R2c = %.3f\n", r2[1], r2[2]))

  aov_tbl <- anova(model)
  aov_tbl$eta2p <- with(aov_tbl, (`F value` * NumDF) / (`F value` * NumDF + DenDF))
  cat("\n-- Type III ANOVA (Satterthwaite) + partial eta^2 --\n")
  print(round(aov_tbl, 4))

  aov_row <- as.data.frame(aov_tbl)
  omnibus_list[[v]] <- data.frame(
    outcome = v, effect = "period",
    F_value = round(aov_row$`F value`, 3),
    df_num = aov_row$NumDF, df_den = round(aov_row$DenDF, 1),
    p_value = round(aov_row$`Pr(>F)`, 4),
    eta2p = round(aov_row$eta2p, 3),
    R2m = round(r2[1], 3), R2c = round(r2[2], 3)
  )

  emm <- emmeans(model, ~ period)
  cat("\n-- Estimated marginal means --\n")
  print(emm)

  vs_bl <- contrast(emm, method = "trt.vs.ctrl", ref = 1, adjust = "bonferroni")
  vs_bl_d <- add_cohens_d(vs_bl, model)
  cat("\n-- Each period vs Baseline (Bonferroni) + Cohen's d --\n")
  print(vs_bl_d)
  vs_bl_d$outcome <- v
  pairwise_list[[v]] <- vs_bl_d

  all_pw <- pairs(emm, adjust = "tukey")
  all_pw_d <- add_cohens_d(all_pw, model)
  cat("\n-- ALL pairwise period comparisons (Tukey) + Cohen's d --\n")
  print(all_pw_d)
  all_pw_d$outcome <- v
  all_pairwise_list[[v]] <- all_pw_d

  df_wide <- df %>% select(participant, period, all_of(v)) %>%
    pivot_wider(names_from = period, values_from = all_of(v))

  bf_tbl <- bind_rows(
    paired_bf(df_wide, "Baseline TMSEEG", "Dosing Session",     "Dosing vs Baseline"),
    paired_bf(df_wide, "Baseline TMSEEG", "1-week follow-up",   "1-Wk vs Baseline"),
    paired_bf(df_wide, "Baseline TMSEEG", "2-week follow-up",   "2-Wk vs Baseline"),
    paired_bf(df_wide, "Baseline TMSEEG", "1-month follow-up",  "1-Mo vs Baseline")
  )
  cat("\n-- Bayes Factors vs Baseline (paired, JZS r = 0.707) --\n")
  print(bf_tbl, row.names = FALSE)
  bf_tbl$outcome <- v
  bf_list[[v]] <- bf_tbl
}

# ── 3. Save summary tables ───────────────────────────────────────────────────
omnibus_summary <- bind_rows(omnibus_list)
write.csv(omnibus_summary,
          "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults/gonogo_mixed_model_omnibus.csv",
          row.names = FALSE)

pairwise_summary <- bind_rows(pairwise_list) %>%
  relocate(outcome) %>%
  rename(p_bonf = p.value) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
write.csv(pairwise_summary,
          "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults/gonogo_vs_baseline_pairwise.csv",
          row.names = FALSE)

all_pairwise_summary <- bind_rows(all_pairwise_list) %>%
  relocate(outcome) %>%
  rename(p_tukey = p.value) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
write.csv(all_pairwise_summary,
          "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults/gonogo_all_pairwise_tukey.csv",
          row.names = FALSE)

bf_summary <- bind_rows(bf_list) %>% relocate(outcome)
write.csv(bf_summary,
          "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults/gonogo_bayes_factors.csv",
          row.names = FALSE)

cat("\n\n", banner, "\nOMNIBUS SUMMARY (period effect, all outcomes)\n", banner, "\n", sep = "")
print(omnibus_summary, row.names = FALSE)

cat("\n\n", banner, "\nSIGNIFICANT PAIRWISE COMPARISONS vs BASELINE (Bonferroni, alpha=.05)\n", banner, "\n", sep = "")
sig <- pairwise_summary %>% filter(p_bonf < 0.05) %>% arrange(outcome, p_bonf)
if (nrow(sig) == 0) cat("None survived Bonferroni correction.\n") else print(sig, row.names = FALSE)

cat("\n\n", banner, "\nSIGNIFICANT ALL-PAIRWISE COMPARISONS, ALL SESSION PAIRS (Tukey, alpha=.05)\n", banner, "\n", sep = "")
sig_all <- all_pairwise_summary %>% filter(p_tukey < 0.05) %>% arrange(outcome, p_tukey)
if (nrow(sig_all) == 0) cat("None survived Tukey correction.\n") else print(sig_all, row.names = FALSE)

cat("\n\n", banner, "\nBAYES FACTORS: STRONG EVIDENCE (BF10 > 3 or BF10 < 1/3)\n", banner, "\n", sep = "")
strong_bf <- bf_summary %>% filter(!is.na(BF10), BF10 > 3 | BF10 < (1/3)) %>% arrange(outcome, desc(BF10))
if (nrow(strong_bf) == 0) cat("None met BF10 > 3 or < 1/3 threshold.\n") else print(strong_bf, row.names = FALSE)

message("\nSaved: gonogo_descriptives.csv, gonogo_mixed_model_omnibus.csv, gonogo_vs_baseline_pairwise.csv, ",
        "gonogo_all_pairwise_tukey.csv, gonogo_bayes_factors.csv")
