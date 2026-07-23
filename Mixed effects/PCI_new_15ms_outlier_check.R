library(tidyverse)
library(lme4)
library(lmerTest)

df <- read.csv("/Users/zsayali1/Documents/OPTE/OPTE/Data/PCI_new_15ms.csv",
               fileEncoding = "UTF-8-BOM") %>%
  mutate(PCI_ST = ifelse(is.nan(PCI_ST), NA, PCI_ST))

keep_sessions <- c("Baseline TMSEEG", "Dosing 1", "Dosing 2", "Dosing 3", "Dosing 4")
timepoint_labels <- c("Baseline TMSEEG"="Baseline","Dosing 1"="90 min",
                      "Dosing 2"="180 min","Dosing 3"="270 min","Dosing 4"="360 min")

df <- df %>%
  filter(Session %in% keep_sessions) %>%
  mutate(
    Timepoint = factor(recode(Session, !!!timepoint_labels),
                       levels = c("Baseline","90 min","180 min","270 min","360 min")),
    Session = factor(recode(Session,
      "Baseline TMSEEG"="Session0","Dosing 1"="Session1","Dosing 2"="Session2",
      "Dosing 3"="Session3","Dosing 4"="Session4"),
      levels = paste0("Session", 0:4)),
    Subject = factor(Subject)
  ) %>%
  filter(!is.na(PCI_ST)) %>%
  group_by(Subject) %>%
  mutate(z_PCI = as.numeric(scale(PCI_ST))) %>%
  ungroup()

cat("\n══════════════════════════════════════════════════════\n")
cat(" Raw PCI_ST values by participant and timepoint\n")
cat("══════════════════════════════════════════════════════\n")
df %>%
  select(Subject, Timepoint, PCI_ST) %>%
  pivot_wider(names_from = Timepoint, values_from = PCI_ST) %>%
  print()

cat("\n══════════════════════════════════════════════════════\n")
cat(" Descriptives per timepoint (raw PCI_ST)\n")
cat("══════════════════════════════════════════════════════\n")
df %>%
  group_by(Timepoint) %>%
  summarise(n=n(), Mean=round(mean(PCI_ST),1), SD=round(sd(PCI_ST),1),
            Min=round(min(PCI_ST),1), Max=round(max(PCI_ST),1), .groups="drop") %>%
  print()

model <- lmer(z_PCI ~ Session + (1 | Subject), data = df, REML = TRUE)
df$resid_std <- scale(residuals(model))[,1]

cat("\n══════════════════════════════════════════════════════\n")
cat(" Standardised residuals > |2.5| (potential outliers)\n")
cat("══════════════════════════════════════════════════════\n")
flagged <- df %>% filter(abs(resid_std) > 2.5) %>%
  select(Subject, Timepoint, PCI_ST, z_PCI, resid_std)
if (nrow(flagged) == 0) cat(" None detected.\n") else print(flagged)
cat("\n Residual summary:\n")
print(summary(df$resid_std))

cat("\n══════════════════════════════════════════════════════\n")
cat(" Leave-one-subject-out: Session F and p\n")
cat("══════════════════════════════════════════════════════\n")
subjects <- levels(df$Subject)
loo_results <- purrr::map_dfr(subjects, function(s) {
  df_loo <- filter(df, Subject != s)
  m <- lmer(z_PCI ~ Session + (1 | Subject), data = df_loo, REML = TRUE)
  a <- anova(m)
  data.frame(Excluded = s,
             F_value  = round(a$`F value`, 3),
             p_value  = round(a$`Pr(>F)`, 4),
             eta2p    = round((a$`F value` * a$NumDF) /
                              (a$`F value` * a$NumDF + a$DenDF), 3))
})
print(loo_results, row.names = FALSE)

p1 <- ggplot(df, aes(x = Timepoint, y = PCI_ST, group = Subject, color = Subject)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2, alpha = 0.8) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               linewidth = 1.5, color = "black", linetype = "dashed") +
  stat_summary(aes(group = 1), fun = mean, geom = "point",
               size = 4, color = "black") +
  labs(title = "PCIst (15ms): Individual Trajectories",
       subtitle = "Colored lines = individual participants; dashed black = group mean",
       x = "Timepoint", y = "PCIst (raw)", color = "Subject") +
  theme_classic(base_size = 13) +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10))

ggsave("/Users/zsayali1/Documents/OPTE/OPTE/Plots/PCI_new_15ms_spaghetti.png",
       p1, width = 8, height = 5, dpi = 150)
message("Saved: Plots/PCI_new_15ms_spaghetti.png")
