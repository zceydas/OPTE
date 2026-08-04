## N-back: group-mean trend plots (all outcome metrics as subplots), with
## n-back level (1/2/3-back) overlaid as separate colored lines on the same
## panel, within-subject SEM error bars (computed separately per level), and
## significance asterisks (vs Baseline, Bonferroni, within each level) pulled
## from nback_vs_baseline_pairwise_by_level.csv.
##
## Data note: participant 16's 1-week / 2-week / 1-month sessions were all-
## zero-response (equipment/attention failure) and are excluded here, per the
## same cleaning rule used in nback_mixed_model.R.

library(tidyverse)

results_dir <- "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults"
plots_dir   <- "/Users/zsayali1/Documents/OPTE/OPTE/Plots"

df <- read.csv(file.path(results_dir, "nback_all_participant_metrics.csv"))
pw <- read.csv(file.path(results_dir, "nback_vs_baseline_pairwise_by_level.csv"))

follow_up_periods <- c("1-week follow-up", "2-week follow-up", "1-month follow-up")
df <- df %>% filter(!(participant == 16 & period %in% follow_up_periods))

period_levels <- c("Baseline TMSEEG", "Dosing Session", "1-week follow-up",
                    "2-week follow-up", "1-month follow-up")
period_labels <- c("Baseline", "Dosing", "1-Wk", "2-Wk", "1-Mo")

outcomes <- c("hit_rate_adj", "fa_rate_adj", "zH", "zFA", "dprime",
              "mean_rt_hit", "mean_rt_fa")

df <- df %>%
  mutate(participant = factor(participant),
         period = factor(period, levels = period_levels, labels = period_labels),
         level  = factor(level, levels = c(1, 2, 3)))

sig_stars <- function(p) case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "")
level_colors <- c("1" = "#2166AC", "2" = "#41AB5D", "3" = "#B2182B")
dodge <- position_dodge(width = 0.3)

# ── Long format ───────────────────────────────────────────────────────────
long_df <- df %>%
  select(participant, period, level, all_of(outcomes)) %>%
  pivot_longer(all_of(outcomes), names_to = "outcome", values_to = "value") %>%
  mutate(outcome = factor(outcome, levels = outcomes))

J <- nlevels(long_df$period)

# Within-subject SEM computed separately WITHIN each outcome x level
p_means <- long_df %>%
  group_by(outcome, level, participant) %>%
  summarise(p_mean = mean(value, na.rm = TRUE), .groups = "drop")
grand_means <- long_df %>%
  group_by(outcome, level) %>%
  summarise(grand_mean = mean(value, na.rm = TRUE), .groups = "drop")

norm_df <- long_df %>%
  left_join(p_means, by = c("outcome", "level", "participant")) %>%
  left_join(grand_means, by = c("outcome", "level")) %>%
  mutate(normalized = value - p_mean + grand_mean)

summary_df <- norm_df %>%
  group_by(outcome, level, period) %>%
  summarise(
    Mean   = mean(value, na.rm = TRUE),
    n      = sum(!is.na(value)),
    WS_SEM = (sd(normalized, na.rm = TRUE) / sqrt(n)) * sqrt(J / (J - 1)),
    .groups = "drop"
  ) %>%
  filter(n > 0)

# ── Significance labels (vs Baseline, within each level) ────────────────────
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
    level = factor(level, levels = c(1, 2, 3)),
    stars = sig_stars(p_bonf)
  ) %>%
  filter(stars != "") %>%
  left_join(summary_df, by = c("outcome", "level", "period")) %>%
  left_join(range_df, by = "outcome") %>%
  mutate(y = Mean + WS_SEM + 0.08 * rng)

# ── Plot ─────────────────────────────────────────────────────────────────────
ggplot(summary_df, aes(x = period, y = Mean, color = level, group = level)) +
  geom_line(linewidth = 1, position = dodge) +
  geom_point(size = 2.6, position = dodge) +
  geom_errorbar(aes(ymin = Mean - WS_SEM, ymax = Mean + WS_SEM),
                width = 0.15, linewidth = 0.6, position = dodge) +
  geom_text(data = sig_df, aes(x = period, y = y, label = stars, color = level),
            position = dodge, size = 4.2, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(name = "N-back level", values = level_colors) +
  facet_wrap(~ outcome, scales = "free_y", ncol = 4) +
  labs(
    x = "Time point", y = NULL,
    title = "N-back: Group Means Across Time Point x Level",
    subtitle = paste0("Mean +/- within-subject SEM (Cousineau-Morey, per level), N = 10 (9 for follow-ups)",
                       " | * p<.05, ** p<.01, *** p<.001 vs Baseline within level (Bonferroni)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 9),
    strip.background = element_blank(),
    strip.text    = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(file.path(plots_dir, "nback_metrics_overview.png"), width = 14, height = 8, dpi = 150)
message("Saved: Plots/nback_metrics_overview.png")
