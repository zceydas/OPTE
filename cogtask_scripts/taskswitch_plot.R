## Task-switching / effort-selection: group-mean trend plots (all outcome
## metrics as subplots), with effort_level overlaid as separate colored lines
## on the same panel, within-subject SEM error bars (computed separately per
## effort_level), and significance asterisks (vs Baseline, Bonferroni,
## within each effort_level) pulled from
## taskswitch_vs_baseline_pairwise_by_effort_level.csv.

library(tidyverse)

results_dir <- "/Users/zsayali1/Documents/OPTE/OPTE/Results/CogTaskResults"
plots_dir   <- "/Users/zsayali1/Documents/OPTE/OPTE/Plots"

df <- read.csv(file.path(results_dir, "effort_selection_metrics.csv"))
pw <- read.csv(file.path(results_dir, "taskswitch_vs_baseline_pairwise_by_effort_level.csv"))

df <- df %>% filter(effort_level != 999)

period_levels <- c("Baseline", "1-week", "1-month")
period_labels <- c("Baseline", "1-Wk", "1-Mo")

outcomes <- c(
  "switch_acc_learning", "repeat_acc_learning", "switch_cost_acc_learning",
  "switch_rt_learning", "repeat_rt_learning", "switch_cost_rt_learning",
  "switch_acc_execution", "repeat_acc_execution", "switch_cost_acc_execution",
  "switch_rt_execution", "repeat_rt_execution", "switch_cost_rt_execution",
  "p_select", "decision_rt"
)

df <- df %>%
  mutate(participant  = factor(participant),
         period       = factor(period, levels = period_levels, labels = period_labels),
         effort_level = factor(effort_level, levels = c(1, 2, 3, 4)))

sig_stars <- function(p) case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "")
effort_colors <- c("1" = "#2166AC", "2" = "#41AB5D", "3" = "#F46D43", "4" = "#B2182B")
dodge <- position_dodge(width = 0.35)

# ── Long format ───────────────────────────────────────────────────────────
long_df <- df %>%
  select(participant, period, effort_level, all_of(outcomes)) %>%
  pivot_longer(all_of(outcomes), names_to = "outcome", values_to = "value") %>%
  mutate(outcome = factor(outcome, levels = outcomes))

J <- nlevels(long_df$period)

# Within-subject SEM computed separately WITHIN each outcome x effort_level
p_means <- long_df %>%
  group_by(outcome, effort_level, participant) %>%
  summarise(p_mean = mean(value, na.rm = TRUE), .groups = "drop")
grand_means <- long_df %>%
  group_by(outcome, effort_level) %>%
  summarise(grand_mean = mean(value, na.rm = TRUE), .groups = "drop")

norm_df <- long_df %>%
  left_join(p_means, by = c("outcome", "effort_level", "participant")) %>%
  left_join(grand_means, by = c("outcome", "effort_level")) %>%
  mutate(normalized = value - p_mean + grand_mean)

summary_df <- norm_df %>%
  group_by(outcome, effort_level, period) %>%
  summarise(
    Mean   = mean(value, na.rm = TRUE),
    n      = sum(!is.na(value)),
    WS_SEM = (sd(normalized, na.rm = TRUE) / sqrt(n)) * sqrt(J / (J - 1)),
    .groups = "drop"
  )

# ── Significance labels (vs Baseline, within each effort_level) ─────────────
range_df <- summary_df %>%
  group_by(outcome) %>%
  summarise(rng = max(Mean + WS_SEM) - min(Mean - WS_SEM), .groups = "drop")

sig_df <- pw %>%
  mutate(
    period = str_remove(contrast, " - Baseline$"),
    period = str_remove_all(period, "^\\(|\\)$"),
    period = recode(period, !!!setNames(period_labels, period_levels)),
    period = factor(period, levels = period_labels),
    outcome = factor(outcome, levels = outcomes),
    effort_level = factor(effort_level, levels = c(1, 2, 3, 4)),
    stars = sig_stars(p_bonf)
  ) %>%
  filter(stars != "") %>%
  left_join(summary_df, by = c("outcome", "effort_level", "period")) %>%
  left_join(range_df, by = "outcome") %>%
  mutate(y = Mean + WS_SEM + 0.08 * rng)

# ── Plot ─────────────────────────────────────────────────────────────────────
ggplot(summary_df, aes(x = period, y = Mean, color = effort_level, group = effort_level)) +
  geom_line(linewidth = 1, position = dodge) +
  geom_point(size = 2.4, position = dodge) +
  geom_errorbar(aes(ymin = Mean - WS_SEM, ymax = Mean + WS_SEM),
                width = 0.15, linewidth = 0.6, position = dodge) +
  geom_text(data = sig_df, aes(x = period, y = y, label = stars, color = effort_level),
            position = dodge, size = 4, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(name = "Effort level", values = effort_colors) +
  facet_wrap(~ outcome, scales = "free_y", ncol = 4) +
  labs(
    x = "Time point", y = NULL,
    title = "Task-Switching / Effort-Selection: Group Means Across Time Point x Effort Level",
    subtitle = paste0("Mean +/- within-subject SEM (Cousineau-Morey, per effort level), N = ",
                       n_distinct(df$participant),
                       " | * p<.05, ** p<.01, *** p<.001 vs Baseline within effort level (Bonferroni)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 9),
    strip.background = element_blank(),
    strip.text    = element_text(face = "bold", size = 8),
    legend.position = "top"
  )

ggsave(file.path(plots_dir, "taskswitch_metrics_overview.png"), width = 15, height = 12, dpi = 150)
message("Saved: Plots/taskswitch_metrics_overview.png")
