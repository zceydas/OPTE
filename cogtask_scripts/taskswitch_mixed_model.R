## Task-switching / effort-selection: descriptive stats, mixed-effects models,
## pairwise comparisons vs. Baseline, and Bayes Factors, for each outcome
## measure across time point and effort level.
##
## Random effect: participant (random intercept)
## Model form: outcome ~ period * effort_level + (1 | participant)
##
## Data note: effort_level == 999 marks catch/no-selection trials (NaN on
## every accuracy/RT/switch-cost column, n_offered == 0). These rows are
## dropped before analysis, per instructions.

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)
library(BayesFactor)

# ── Load & prepare ───────────────────────────────────────────────────────────
df <- read.csv("/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults/effort_selection_metrics.csv")

df <- df %>% filter(effort_level != 999)

period_levels <- c("Baseline", "1-week", "1-month")

df <- df %>%
  mutate(
    participant  = factor(participant),
    period       = factor(period, levels = period_levels),
    effort_level = factor(effort_level, levels = c(1, 2, 3, 4))
  )

stopifnot(sum(is.na(df)) == 0)

outcomes <- c(
  "switch_acc_learning", "repeat_acc_learning", "switch_cost_acc_learning",
  "switch_rt_learning", "repeat_rt_learning", "switch_cost_rt_learning",
  "switch_acc_execution", "repeat_acc_execution", "switch_cost_acc_execution",
  "switch_rt_execution", "repeat_rt_execution", "switch_cost_rt_execution",
  "p_select", "decision_rt"
)

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
# Marginal (collapsed across effort_level)
desc_marg_list <- list()
for (v in outcomes) {
  desc_marg_list[[v]] <- df %>%
    group_by(period) %>%
    summarise(outcome = v, n = sum(!is.na(.data[[v]])),
              mean = mean(.data[[v]], na.rm = TRUE),
              sd = sd(.data[[v]], na.rm = TRUE), .groups = "drop") %>%
    relocate(outcome)
}
descriptives_marginal <- bind_rows(desc_marg_list) %>%
  mutate(across(c(mean, sd), ~ round(.x, 4)))

write.csv(descriptives_marginal,
          "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults/taskswitch_descriptives_marginal.csv",
          row.names = FALSE)

# Cell-level (period x effort_level)
desc_cell_list <- list()
for (v in outcomes) {
  desc_cell_list[[v]] <- df %>%
    group_by(period, effort_level) %>%
    summarise(outcome = v, n = sum(!is.na(.data[[v]])),
              mean = mean(.data[[v]], na.rm = TRUE),
              sd = sd(.data[[v]], na.rm = TRUE), .groups = "drop") %>%
    relocate(outcome)
}
descriptives_cell <- bind_rows(desc_cell_list) %>%
  mutate(across(c(mean, sd), ~ round(.x, 4)))

write.csv(descriptives_cell,
          "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults/taskswitch_descriptives_by_effort_level.csv",
          row.names = FALSE)

cat("\n", banner, "\nDESCRIPTIVE STATISTICS -- marginal (collapsed over effort_level)\n", banner, "\n", sep = "")
print(descriptives_marginal, n = Inf)

# ── 2. Mixed models, omnibus, pairwise vs Baseline, Bayes Factors ───────────
omnibus_list          <- list()
pairwise_marg_list    <- list()
pairwise_by_eff_list  <- list()
bf_list               <- list()

for (v in outcomes) {

  cat("\n\n", banner, "\nOUTCOME: ", v, "\n", banner, "\n", sep = "")

  form  <- as.formula(paste0(v, " ~ period * effort_level + (1 | participant)"))
  model <- lmer(form, data = df, REML = TRUE, na.action = na.omit)

  r2 <- r.squaredGLMM(model)
  cat(sprintf("\nR2m = %.3f | R2c = %.3f\n", r2[1], r2[2]))

  aov_tbl <- anova(model)
  aov_tbl$eta2p <- with(aov_tbl, (`F value` * NumDF) / (`F value` * NumDF + DenDF))
  cat("\n-- Type III ANOVA (Satterthwaite) + partial eta^2 --\n")
  print(round(aov_tbl, 4))

  aov_df <- as.data.frame(aov_tbl)
  aov_df$effect <- rownames(aov_df)
  omnibus_list[[v]] <- data.frame(
    outcome = v, effect = aov_df$effect,
    F_value = round(aov_df$`F value`, 3),
    df_num = aov_df$NumDF, df_den = round(aov_df$DenDF, 1),
    p_value = round(aov_df$`Pr(>F)`, 4),
    eta2p = round(aov_df$eta2p, 3),
    R2m = round(r2[1], 3), R2c = round(r2[2], 3)
  )

  # -- Marginal period comparisons (averaged over effort_level) --
  emm_marg <- emmeans(model, ~ period)
  vs_bl_marg <- contrast(emm_marg, method = "trt.vs.ctrl", ref = 1, adjust = "bonferroni")
  vs_bl_marg_d <- add_cohens_d(vs_bl_marg, model)
  cat("\n-- Period vs Baseline, averaged over effort_level (Bonferroni) + Cohen's d --\n")
  print(vs_bl_marg_d)
  vs_bl_marg_d$outcome <- v
  pairwise_marg_list[[v]] <- vs_bl_marg_d

  # -- Period comparisons within each effort_level --
  emm_eff <- emmeans(model, ~ period | effort_level)
  vs_bl_eff <- contrast(emm_eff, method = "trt.vs.ctrl", ref = 1, adjust = "bonferroni")
  vs_bl_eff_d <- add_cohens_d(vs_bl_eff, model)
  cat("\n-- Period vs Baseline, within each effort_level (Bonferroni) + Cohen's d --\n")
  print(vs_bl_eff_d)
  vs_bl_eff_d$outcome <- v
  pairwise_by_eff_list[[v]] <- vs_bl_eff_d

  # -- Bayes Factors vs Baseline (subject-level means, collapsed over effort_level) --
  df_wide <- df %>%
    group_by(participant, period) %>%
    summarise(m = mean(.data[[v]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = period, values_from = m)

  bf_tbl <- bind_rows(
    paired_bf(df_wide, "Baseline", "1-week",  "1-Wk vs Baseline"),
    paired_bf(df_wide, "Baseline", "1-month", "1-Mo vs Baseline")
  )
  cat("\n-- Bayes Factors vs Baseline (paired, subject means, JZS r = 0.707) --\n")
  print(bf_tbl, row.names = FALSE)
  bf_tbl$outcome <- v
  bf_list[[v]] <- bf_tbl
}

# ── 3. Save summary tables ───────────────────────────────────────────────────
omnibus_summary <- bind_rows(omnibus_list)
write.csv(omnibus_summary,
          "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults/taskswitch_mixed_model_omnibus.csv",
          row.names = FALSE)

pairwise_marg_summary <- bind_rows(pairwise_marg_list) %>%
  relocate(outcome) %>% rename(p_bonf = p.value) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
write.csv(pairwise_marg_summary,
          "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults/taskswitch_vs_baseline_pairwise_marginal.csv",
          row.names = FALSE)

pairwise_by_eff_summary <- bind_rows(pairwise_by_eff_list) %>%
  relocate(outcome) %>% rename(p_bonf = p.value) %>%
  mutate(across(where(is.numeric) & !any_of("effort_level"), ~ round(.x, 4)))
write.csv(pairwise_by_eff_summary,
          "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults/taskswitch_vs_baseline_pairwise_by_effort_level.csv",
          row.names = FALSE)

bf_summary <- bind_rows(bf_list) %>% relocate(outcome)
write.csv(bf_summary,
          "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults/taskswitch_bayes_factors.csv",
          row.names = FALSE)

cat("\n\n", banner, "\nOMNIBUS SUMMARY (period, effort_level, interaction)\n", banner, "\n", sep = "")
print(omnibus_summary, row.names = FALSE)

cat("\n\n", banner, "\nSIGNIFICANT PAIRWISE vs BASELINE, MARGINAL (Bonferroni, alpha=.05)\n", banner, "\n", sep = "")
sig_marg <- pairwise_marg_summary %>% filter(p_bonf < 0.05) %>% arrange(outcome, p_bonf)
if (nrow(sig_marg) == 0) cat("None survived Bonferroni correction.\n") else print(sig_marg, row.names = FALSE)

cat("\n\n", banner, "\nBAYES FACTORS: STRONG EVIDENCE (BF10 > 3 or BF10 < 1/3)\n", banner, "\n", sep = "")
strong_bf <- bf_summary %>% filter(!is.na(BF10), BF10 > 3 | BF10 < (1/3)) %>% arrange(outcome, desc(BF10))
if (nrow(strong_bf) == 0) cat("None met BF10 > 3 or < 1/3 threshold.\n") else print(strong_bf, row.names = FALSE)

message("\nSaved: taskswitch_descriptives_marginal.csv, taskswitch_descriptives_by_effort_level.csv, ",
        "taskswitch_mixed_model_omnibus.csv, taskswitch_vs_baseline_pairwise_marginal.csv, ",
        "taskswitch_vs_baseline_pairwise_by_effort_level.csv, taskswitch_bayes_factors.csv")
