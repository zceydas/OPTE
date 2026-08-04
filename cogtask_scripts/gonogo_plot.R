## Go/No-Go: group-mean trend plots (all outcome metrics as subplots), with
## within-subject SEM error bars and significance asterisks (vs Baseline,
## Bonferroni-corrected) pulled from gonogo_vs_baseline_pairwise.csv.

library(tidyverse)

results_dir <- "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults"
plots_dir   <- "/Users/zsayali1/Documents/OPTE/OPTE/Plots"

df <- read.csv(file.path(results_dir, "gonogo_metrics_by_period.csv"))
pw <- read.csv(file.path(results_dir, "gonogo_vs_baseline_pairwise.csv"))

period_levels <- c("Baseline TMSEEG", "Dosing Session", "1-week follow-up",
                    "2-week follow-up", "1-month follow-up")
period_labels <- c("Baseline", "Dosing", "1-Wk", "2-Wk", "1-Mo")

outcomes <- c("H", "FA", "zH", "zFA", "dprime", "criterion", "mean_rt")

df <- df %>%
  mutate(participant = factor(participant),
         period = factor(period, levels = period_levels, labels = period_labels))

sig_stars <- function(p) case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "")

# ── Long format + within-subject (Cousineau-Morey) SEM per outcome ──────────
long_df <- df %>%
  select(participant, period, all_of(outcomes)) %>%
  pivot_longer(all_of(outcomes), names_to = "outcome", values_to = "value") %>%
  mutate(outcome = factor(outcome, levels = outcomes))

J <- nlevels(long_df$period)

p_means <- long_df %>%
  group_by(outcome, participant) %>%
  summarise(p_mean = mean(value, na.rm = TRUE), .groups = "drop")
grand_means <- long_df %>%
  group_by(outcome) %>%
  summarise(grand_mean = mean(value, na.rm = TRUE), .groups = "drop")

norm_df <- long_df %>%
  left_join(p_means, by = c("outcome", "participant")) %>%
  left_join(grand_means, by = "outcome") %>%
  mutate(normalized = value - p_mean + grand_mean)

summary_df <- norm_df %>%
  group_by(outcome, period) %>%
  summarise(
    Mean   = mean(value, na.rm = TRUE),
    n      = sum(!is.na(value)),
    WS_SEM = (sd(normalized, na.rm = TRUE) / sqrt(n)) * sqrt(J / (J - 1)),
    .groups = "drop"
  )

# ── Significance labels (vs Baseline) ────────────────────────────────────────
range_df <- summary_df %>%
  group_by(outcome) %>%
  summarise(rng = max(Mean + WS_SEM) - min(Mean - WS_SEM), .groups = "drop")

sig_df <- pw %>%
  mutate(
    period = str_remove(contrast, " - Baseline TMSEEG$"),
    period = str_remove_all(period, "^\\(|\\)$"),
    period = recode(period, !!!setNames(period_labels, period_levels)),
    period = factor(period, levels = period_labels),
    outcome = factor(outcome, levels = outcomes),
    stars = sig_stars(p_bonf)
  ) %>%
  filter(stars != "") %>%
  left_join(summary_df, by = c("outcome", "period")) %>%
  left_join(range_df, by = "outcome") %>%
  mutate(y = Mean + WS_SEM + 0.08 * rng)

# ── Plot ─────────────────────────────────────────────────────────────────────
ggplot(summary_df, aes(x = period, y = Mean, group = 1)) +
  geom_line(linewidth = 1.1, color = "#2166AC") +
  geom_point(size = 2.8, color = "#2166AC") +
  geom_errorbar(aes(ymin = Mean - WS_SEM, ymax = Mean + WS_SEM),
                width = 0.15, linewidth = 0.7, color = "#2166AC") +
  geom_text(data = sig_df, aes(x = period, y = y, label = stars),
            inherit.aes = FALSE, size = 5, color = "black", fontface = "bold") +
  facet_wrap(~ outcome, scales = "free_y", ncol = 4) +
  labs(
    x = "Time point", y = NULL,
    title = "Go/No-Go: Group Means Across Time Point",
    subtitle = paste0("Mean +/- within-subject SEM (Cousineau-Morey), N = ",
                       n_distinct(df$participant),
                       " | * p<.05, ** p<.01, *** p<.001 vs Baseline (Bonferroni)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 9),
    strip.background = element_blank(),
    strip.text    = element_text(face = "bold")
  )

ggsave(file.path(plots_dir, "gonogo_metrics_overview.png"), width = 13, height = 7, dpi = 150)
message("Saved: Plots/gonogo_metrics_overview.png")
