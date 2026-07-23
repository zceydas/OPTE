library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(effectsize)
library(MuMIn)
library(BayesFactor)

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
    Session = factor(Session, levels = paste0("Session", 0:4)),
    Subject = factor(Subject)
  ) %>%
  filter(!is.na(PCI_ST))

df <- df %>%
  group_by(Subject) %>%
  mutate(z_PCI = as.numeric(scale(PCI_ST))) %>%
  ungroup()

best_model <- lmer(z_PCI ~ Session + (1 | Subject), data = df, REML = TRUE)

r2_vals <- r.squaredGLMM(best_model)
cat("\n══════════════════════════════════════════════════════\n")
cat(" R² (marginal / conditional)\n")
cat("══════════════════════════════════════════════════════\n")
print(round(r2_vals, 3))

aov_tbl <- anova(best_model)
aov_tbl$eta2p <- with(aov_tbl, (`F value` * NumDF) / (`F value` * NumDF + DenDF))
cat("\n══════════════════════════════════════════════════════\n")
cat(" Type III ANOVA + partial η²\n")
cat("══════════════════════════════════════════════════════\n")
print(round(aov_tbl, 4))

cat("\n══════════════════════════════════════════════════════\n")
cat(" Model summary\n")
cat("══════════════════════════════════════════════════════\n")
print(summary(best_model))

emm <- emmeans(best_model, ~ Session)
cat("\n══════════════════════════════════════════════════════\n")
cat(" Estimated marginal means (z_PCI)\n")
cat("══════════════════════════════════════════════════════\n")
print(emm)

add_cohens_d <- function(contrast_obj, model) {
  tbl <- as.data.frame(contrast_obj)
  tbl$d <- tbl$estimate / sigma(model)
  tbl
}

session_contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = 1,
                               adjust = "fdr")
cat("\n══════════════════════════════════════════════════════\n")
cat(" Each Timepoint vs Baseline (FDR) + Cohen's d\n")
cat("══════════════════════════════════════════════════════\n")
print(add_cohens_d(session_contrasts, best_model))

pairwise_contrasts <- contrast(emm, method = "pairwise", adjust = "fdr")
cat("\n══════════════════════════════════════════════════════\n")
cat(" All pairwise comparisons (FDR) + Cohen's d\n")
cat("══════════════════════════════════════════════════════\n")
print(add_cohens_d(pairwise_contrasts, best_model))

df_wide <- df %>%
  select(Subject, Session, PCI_ST) %>%
  pivot_wider(names_from = Session, values_from = PCI_ST)

bf_table <- function(wide, s1, s2, label) {
  d <- wide[[s2]] - wide[[s1]]
  d <- d[!is.na(d)]
  n  <- length(d)
  dz <- mean(d) / sd(d)
  bf <- ttestBF(x = d, mu = 0)
  bf10 <- exp(bf@bayesFactor$bf)
  data.frame(Contrast=label, n_pairs=n, mean_diff=round(mean(d),3),
             dz=round(dz,3), BF10=round(bf10,3), BF01=round(1/bf10,3))
}

vs_baseline <- bind_rows(
  bf_table(df_wide, "Session0", "Session1", "90 min vs Baseline"),
  bf_table(df_wide, "Session0", "Session2", "180 min vs Baseline"),
  bf_table(df_wide, "Session0", "Session3", "270 min vs Baseline"),
  bf_table(df_wide, "Session0", "Session4", "360 min vs Baseline")
)
cat("\n══════════════════════════════════════════════════════\n")
cat(" Bayes Factors: Each Timepoint vs Baseline (JZS, r = 0.707)\n")
cat("══════════════════════════════════════════════════════\n")
print(vs_baseline, row.names = FALSE)

df <- df %>% mutate(session_num = as.numeric(Session) - 1)

poly_model <- lmer(z_PCI ~ poly(session_num, 2) + (1 | Subject),
                   data = df, REML = TRUE)
aov_poly <- anova(poly_model)
aov_poly$eta2p <- with(aov_poly, (`F value` * NumDF) / (`F value` * NumDF + DenDF))
r2_poly <- r.squaredGLMM(poly_model)

cat("\n══════════════════════════════════════════════════════\n")
cat(" Polynomial Trend Model: z_PCI ~ poly(Session, 2)\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("\n R²m = %.3f  |  R²c = %.3f\n\n", r2_poly[1], r2_poly[2]))
print(round(aov_poly, 4))
cat("\n Fixed-effect coefficients\n")
print(round(as.data.frame(summary(poly_model)$coefficients), 4))
