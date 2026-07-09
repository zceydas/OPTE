library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(effectsize)  # eta_squared, interpret_eta_squared
library(MuMIn)       # r.squaredGLMM

# ── Load & prepare ─────────────────────────────────────────────────────────────
df <- read.csv("/Users/zsayali1/Documents/OPTE/OPTE/Data/GLOBAL_MEAN_LZc_SUMMARY.csv")

# Within-subject z-score. Rationale: between-subject SD (0.063) is nearly as
# large as within-subject SD (0.088); z-scoring removes individual differences
# in absolute LZc level and scale, leaving only the within-person structure.
df <- df %>%
  group_by(Participant) %>%
  mutate(z_LZc = as.numeric(scale(Global_Mean_LZc))) %>%
  ungroup() %>%
  mutate(
    Participant = factor(Participant),
    Eyes        = factor(Eyes,  levels = c("EC", "EO")),
    Epoch       = factor(Epoch, levels = paste0("Epoch", 0:4))
  )

# ── Fit model ─────────────────────────────────────────────────────────────────
# Random slope for Eyes tested but singular (LRT p = 1); use random intercept.
best_model <- lmer(z_LZc ~ Eyes * Epoch + (1 | Participant),
                   data = df, REML = TRUE)

# ── R² (marginal = fixed effects; conditional = fixed + random) ───────────────
r2_vals <- r.squaredGLMM(best_model)

cat("\n══════════════════════════════════════════════════════\n")
cat(" R² (marginal / conditional)\n")
cat("══════════════════════════════════════════════════════\n")
print(round(r2_vals, 3))

# ── Type III ANOVA + partial η² ───────────────────────────────────────────────
aov_tbl <- anova(best_model)   # lmerTest: Satterthwaite F + p

# partial η²p = (F * df_num) / (F * df_num + df_den)  [computed from F & df]
aov_tbl$eta2p <- with(aov_tbl,
  (`F value` * NumDF) / (`F value` * NumDF + DenDF)
)

cat("\n══════════════════════════════════════════════════════\n")
cat(" Type III ANOVA + partial η²\n")
cat("══════════════════════════════════════════════════════\n")
print(round(aov_tbl, 4))

# ── Model summary ─────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════\n")
cat(" Model summary\n")
cat("══════════════════════════════════════════════════════\n")
print(summary(best_model))

# ── Marginal means ────────────────────────────────────────────────────────────
emm <- emmeans(best_model, ~ Eyes * Epoch)

cat("\n══════════════════════════════════════════════════════\n")
cat(" Estimated marginal means (z_LZc)\n")
cat("══════════════════════════════════════════════════════\n")
print(emm)

# Cohen's d helper: d = estimate / residual sigma (within-person SD after
# removing fixed effects). Appropriate here because data are z-scored and
# sigma captures the within-participant residual variability.
add_cohens_d <- function(contrast_obj, model) {
  s <- sigma(model)
  tbl <- as.data.frame(contrast_obj)
  tbl$d <- tbl$estimate / s
  tbl
}

# ── 1. Main effect of Eyes (collapsed across epochs) ─────────────────────────
emm_eyes  <- emmeans(best_model, ~ Eyes)
eyes_main <- contrast(emm_eyes, method = "pairwise")

cat("\n══════════════════════════════════════════════════════\n")
cat(" 1. Main effect of Eyes (averaged across epochs) + Cohen's d\n")
cat("══════════════════════════════════════════════════════\n")
print(add_cohens_d(eyes_main, best_model))

# ── 2. EC vs EO at each Epoch (Bonferroni corrected for 5 tests) ──────────────
# `by = "Epoch"` runs 1 test per epoch; manually apply Bonferroni × 5 across
# all five EC–EO comparisons so they form a single corrected family.
eyes_by_epoch_raw <- as.data.frame(
  contrast(emm, method = "pairwise", by = "Epoch")
)
eyes_by_epoch_raw$p_bonf <- pmin(eyes_by_epoch_raw$p.value * 5, 1)
eyes_by_epoch_raw$d      <- eyes_by_epoch_raw$estimate / sigma(best_model)

cat("\n══════════════════════════════════════════════════════\n")
cat(" 2. EC vs EO at each Epoch (Bonferroni × 5) + Cohen's d\n")
cat("══════════════════════════════════════════════════════\n")
print(eyes_by_epoch_raw[, c("Epoch", "contrast", "estimate", "SE",
                             "df", "t.ratio", "p.value", "p_bonf", "d")])

# ── 3. Each Epoch vs Baseline — separately for EC and EO ─────────────────────
epoch_vs_base <- contrast(emm, method = "trt.vs.ctrl", ref = 1,
                           by = "Eyes", adjust = "bonferroni")

cat("\n══════════════════════════════════════════════════════\n")
cat(" 3. Each Epoch vs Baseline (Epoch0), by condition (Bonferroni × 4) + d\n")
cat("══════════════════════════════════════════════════════\n")
print(add_cohens_d(epoch_vs_base, best_model))

# ── 4. Epoch vs Baseline — collapsed across Eyes ──────────────────────────────
emm_epoch        <- emmeans(best_model, ~ Epoch)
epoch_vs_base_av <- contrast(emm_epoch, method = "trt.vs.ctrl", ref = 1,
                              adjust = "bonferroni")

cat("\n══════════════════════════════════════════════════════\n")
cat(" 4. Each Epoch vs Baseline (averaged across Eyes, Bonferroni × 4) + d\n")
cat("══════════════════════════════════════════════════════\n")
print(add_cohens_d(epoch_vs_base_av, best_model))

# ── 5. All pairwise Epoch comparisons (averaged across Eyes) ──────────────────
# C(5,2) = 10 pairwise comparisons, Bonferroni corrected for 10 tests.
epoch_pairwise <- contrast(emm_epoch, method = "pairwise", adjust = "bonferroni")

cat("\n══════════════════════════════════════════════════════\n")
cat(" 5. All pairwise Epoch comparisons (Bonferroni × 10) + Cohen's d\n")
cat("══════════════════════════════════════════════════════\n")
print(add_cohens_d(epoch_pairwise, best_model))

# ── 6. All pairwise Epoch comparisons — separately for EC and EO ──────────────
epoch_pairwise_by_eyes <- contrast(emm, method = "pairwise",
                                    by = "Eyes", adjust = "bonferroni")

cat("\n══════════════════════════════════════════════════════\n")
cat(" 6. All pairwise Epoch comparisons, by condition (Bonferroni × 10) + d\n")
cat("══════════════════════════════════════════════════════\n")
print(add_cohens_d(epoch_pairwise_by_eyes, best_model))

# ── Polynomial trend model ────────────────────────────────────────────────────
# Fits linear and quadratic time trends (paralleling the MEQ4 analysis) so
# that the shape of the LZc trajectory can be directly compared with PCI and
# MEQ4. Epoch_num = 0–4 (numeric); poly() produces orthogonal polynomials so
# linear and quadratic terms are independent. Eyes interaction tested to confirm
# the trajectory is consistent across EC/EO conditions.

df <- df %>% mutate(Epoch_num = as.numeric(Epoch) - 1)  # 0,1,2,3,4

poly_model <- lmer(z_LZc ~ Eyes * poly(Epoch_num, 2) + (1 | Participant),
                   data = df, REML = TRUE)

aov_poly <- anova(poly_model)
aov_poly$eta2p <- with(aov_poly,
  (`F value` * NumDF) / (`F value` * NumDF + DenDF))
r2_poly <- r.squaredGLMM(poly_model)

cat("\n══════════════════════════════════════════════════════\n")
cat(" LZc — Polynomial Trend Model: z_LZc ~ Eyes × poly(Epoch, 2)\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("\n R²m = %.3f  |  R²c = %.3f\n\n", r2_poly[1], r2_poly[2]))
cat(" Type III ANOVA + partial η²\n")
print(round(aov_poly, 4))

cat("\n Fixed-effect coefficients (linear vs quadratic decomposition)\n")
poly_coefs <- as.data.frame(summary(poly_model)$coefficients)
print(round(poly_coefs, 4))

cat("\n Interpretation:\n")
cat("  poly1 = linear trend across epochs\n")
cat("  poly2 = quadratic (curvature: negative = inverted-U = rise then fall)\n")
cat("  Eyes:poly interaction tests whether trajectory shape differs EC vs EO\n")
