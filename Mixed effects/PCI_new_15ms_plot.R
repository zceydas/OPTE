library(tidyverse)

df <- read.csv("/Users/zsayali1/Documents/OPTE/OPTE/Data/PCI_new_15ms.csv",
               fileEncoding = "UTF-8-BOM") %>%
  mutate(PCI_ST = ifelse(is.nan(PCI_ST), NA, PCI_ST))

keep_sessions <- c("Baseline TMSEEG", "Dosing 1", "Dosing 2", "Dosing 3", "Dosing 4")

df <- df %>%
  filter(Session %in% keep_sessions) %>%
  mutate(
    Session = recode(Session,
      "Baseline TMSEEG" = "Session0",
      "Dosing 1"        = "Session1",
      "Dosing 2"        = "Session2",
      "Dosing 3"        = "Session3",
      "Dosing 4"        = "Session4"
    ),
    Session = factor(Session, levels = paste0("Session", 0:4))
  ) %>%
  filter(!is.na(PCI_ST))

session_labels <- c("Baseline", "90 min", "180 min", "270 min", "360 min")

p_means   <- df %>% group_by(Subject) %>%
  summarise(p_mean = mean(PCI_ST, na.rm = TRUE), .groups = "drop")
grand_mean <- mean(df$PCI_ST, na.rm = TRUE)
J          <- nrow(distinct(df, Session))

df <- df %>%
  left_join(p_means, by = "Subject") %>%
  mutate(normalized = PCI_ST - p_mean + grand_mean)

summary_df <- df %>%
  group_by(Session) %>%
  summarise(
    Mean   = mean(PCI_ST, na.rm = TRUE),
    n      = sum(!is.na(PCI_ST)),
    WS_SEM = (sd(normalized, na.rm = TRUE) / sqrt(n)) * sqrt(J / (J - 1)),
    .groups = "drop"
  )

ggplot(summary_df, aes(x = Session, y = Mean, group = 1)) +
  geom_line(linewidth = 1.2, color = "#2166AC") +
  geom_point(size = 3.5, color = "#2166AC") +
  geom_errorbar(
    aes(ymin = Mean - WS_SEM, ymax = Mean + WS_SEM),
    width = 0.18, linewidth = 0.9, color = "#2166AC"
  ) +
  scale_x_discrete(labels = session_labels) +
  labs(
    x        = "Timepoint",
    y        = "PCIst",
    title    = "PCIst Across Sessions (15ms blanking)",
    subtitle = paste0("Mean ± within-subject SEM (Cousineau–Morey correction), N = ",
                      length(unique(df$Subject)))
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10),
    axis.text.x   = element_text(size = 10)
  )

ggsave("/Users/zsayali1/Documents/OPTE/OPTE/Plots/PCI_new_15ms_across_sessions.png",
       width = 7, height = 5, dpi = 150)
message("Saved: Plots/PCI_new_15ms_across_sessions.png")
