library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)

# ══════════════════════════════════════════════════════════════════════════════
#  Research question: Does within-session LZc predict concurrent MEQ4 scores?
#
#  Timepoint mapping (confirmed by P008 missing-data alignment):
#    MEQ T4 ↔ EEG Epoch1
#    MEQ T6 ↔ EEG Epoch2
#    MEQ T8 ↔ EEG Epoch3
#    MEQ T9 ↔ EEG Epoch4
#
#  LZc: averaged across EC and EO (Eyes × Epoch interaction was non-significant),
#  then z-scored within participant across ALL epochs (0–4) so the scale is
#  consistent with the main LZc analysis and baseline is part of the reference.
# ══════════════════════════════════════════════════════════════════════════════

# ── 1. MEQ4 DATA ──────────────────────────────────────────────────────────────
# Values from dosing session; T4, T6, T8, T9 are the concurrent timepoints
meq_data <- tribble(
  ~Participant, ~T4, ~T6, ~T8, ~T9,
   5,            1,   3,   1,   1,
   7,            3,   3,   4,   2,
   8,           NA,   5,   2,   2,
  11,            0,   1,   1,   0,
  12,            5,   5,   3,   3,
  13,            3,   5,   0,   0,
  16,            3,   3,   1,   1,
  17,            5,   5,   5,   4,
  18,            2,   4,   4,   4,
  19,            5,   3,   0,   0
)

meq_long <- meq_data %>%
  pivot_longer(cols = c(T4, T6, T8, T9),
               names_to  = "MEQ_timepoint",
               values_to = "MEQ4") %>%
  mutate(
    Epoch    = recode(MEQ_timepoint,
                      T4 = "Epoch1", T6 = "Epoch2",
                      T8 = "Epoch3", T9 = "Epoch4"),
    Epoch    = factor(Epoch, levels = paste0("Epoch", 1:4))
  )

# ── 2. LZc DATA ───────────────────────────────────────────────────────────────
lzc_raw <- read.csv("/Users/zsayali1/Documents/OPTE/OPTE/Data/GLOBAL_MEAN_LZc_SUMMARY.csv")

# Average EC and EO per participant × epoch (interaction was non-significant)
lzc_avg <- lzc_raw %>%
  group_by(Participant, Epoch) %>%
  summarise(Mean_LZc = mean(Global_Mean_LZc, na.rm = TRUE), .groups = "drop") %>%
  mutate(Epoch = factor(Epoch, levels = paste0("Epoch", 0:4)))

# Within-subject z-score across ALL epochs (including baseline) so the
# reference distribution includes the pre-drug state
lzc_z <- lzc_avg %>%
  group_by(Participant) %>%
  mutate(z_LZc = as.numeric(scale(Mean_LZc))) %>%
  ungroup() %>%
  filter(Epoch != "Epoch0")   # keep only dosing epochs for the merge

# ── 3. MERGE ──────────────────────────────────────────────────────────────────
df <- meq_long %>%
  left_join(lzc_z, by = c("Participant", "Epoch")) %>%
  filter(!is.na(MEQ4), !is.na(z_LZc)) %>%
  mutate(Participant = factor(Participant))

cat("\n══════════════════════════════════════════════════════\n")
cat(" Merged dataset\n")
cat("══════════════════════════════════════════════════════\n")
print(df %>% select(Participant, Epoch, MEQ_timepoint, MEQ4, Mean_LZc, z_LZc),
      n = Inf)
cat(sprintf("\nN obs = %d across %d participants\n", nrow(df), n_distinct(df$Participant)))

# ── 4. Cohen's d helper ───────────────────────────────────────────────────────
add_cohens_d <- function(contrast_obj, model) {
  tbl <- as.data.frame(contrast_obj)
  tbl$d <- tbl$estimate / sigma(model)
  tbl
}

# ══════════════════════════════════════════════════════════════════════════════
#  MODEL 1 — Simple association: MEQ4 ~ z_LZc
#  Does LZc predict MEQ4 across the dosing session (no time covariate)?
# ══════════════════════════════════════════════════════════════════════════════
m1 <- lmer(MEQ4 ~ z_LZc + (1 | Participant), data = df, REML = TRUE)

cat("\n══════════════════════════════════════════════════════\n")
cat(" Model 1: MEQ4 ~ z_LZc + (1 | Participant)\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("\n R²m = %.3f  |  R²c = %.3f\n", r.squaredGLMM(m1)[1], r.squaredGLMM(m1)[2]))
print(summary(m1))
aov1 <- anova(m1)
aov1$eta2p <- with(aov1, (`F value` * NumDF) / (`F value` * NumDF + DenDF))
cat(" ANOVA:\n"); print(round(aov1, 4))

# ══════════════════════════════════════════════════════════════════════════════
#  MODEL 2 — Time-controlled: MEQ4 ~ z_LZc + Epoch
#  Does LZc predict MEQ4 above and beyond the shared temporal trend?
#  (Both LZc and MEQ4 show inverted-U trajectories; this controls for that.)
# ══════════════════════════════════════════════════════════════════════════════
m2 <- lmer(MEQ4 ~ z_LZc + Epoch + (1 | Participant), data = df, REML = TRUE)

cat("\n══════════════════════════════════════════════════════\n")
cat(" Model 2: MEQ4 ~ z_LZc + Epoch + (1 | Participant)\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("\n R²m = %.3f  |  R²c = %.3f\n", r.squaredGLMM(m2)[1], r.squaredGLMM(m2)[2]))
print(summary(m2))
aov2 <- anova(m2)
aov2$eta2p <- with(aov2, (`F value` * NumDF) / (`F value` * NumDF + DenDF))
cat(" ANOVA:\n"); print(round(aov2, 4))

# ══════════════════════════════════════════════════════════════════════════════
#  MODEL 3 — Interaction: MEQ4 ~ z_LZc * Epoch
#  Does the LZc→MEQ4 relationship vary across time points?
# ══════════════════════════════════════════════════════════════════════════════
m3 <- lmer(MEQ4 ~ z_LZc * Epoch + (1 | Participant), data = df, REML = TRUE)

cat("\n══════════════════════════════════════════════════════\n")
cat(" Model 3: MEQ4 ~ z_LZc * Epoch + (1 | Participant)\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("\n R²m = %.3f  |  R²c = %.3f\n", r.squaredGLMM(m3)[1], r.squaredGLMM(m3)[2]))
print(summary(m3))
aov3 <- anova(m3)
aov3$eta2p <- with(aov3, (`F value` * NumDF) / (`F value` * NumDF + DenDF))
cat(" ANOVA:\n"); print(round(aov3, 4))

# ── Model comparison (ML for LRT) ─────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════\n")
cat(" Model comparison (likelihood ratio tests, ML)\n")
cat("══════════════════════════════════════════════════════\n")
m1_ml <- lmer(MEQ4 ~ z_LZc             + (1 | Participant), data = df, REML = FALSE)
m2_ml <- lmer(MEQ4 ~ z_LZc + Epoch     + (1 | Participant), data = df, REML = FALSE)
m3_ml <- lmer(MEQ4 ~ z_LZc * Epoch     + (1 | Participant), data = df, REML = FALSE)
print(anova(m1_ml, m2_ml, m3_ml))

# ── Best model: report emmeans & simple slopes ────────────────────────────────
# Use Model 2 (time-controlled) as the primary interpretive model
cat("\n══════════════════════════════════════════════════════\n")
cat(" Model 2 — estimated marginal means of MEQ4 by Epoch\n")
cat("══════════════════════════════════════════════════════\n")
emm2 <- emmeans(m2, ~ Epoch)
print(emm2)

cat("\n══════════════════════════════════════════════════════\n")
cat(" Model 2 — z_LZc effect (Cohen's d)\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("  b = %.4f, SE = %.4f, t = %.3f, p = %.4f\n",
            fixef(m2)["z_LZc"],
            sqrt(vcov(m2)["z_LZc","z_LZc"]),
            summary(m2)$coefficients["z_LZc","t value"],
            summary(m2)$coefficients["z_LZc","Pr(>|t|)"]))
cat(sprintf("  Cohen's d = %.3f  (b / sigma = %.4f / %.4f)\n",
            fixef(m2)["z_LZc"] / sigma(m2),
            fixef(m2)["z_LZc"], sigma(m2)))
