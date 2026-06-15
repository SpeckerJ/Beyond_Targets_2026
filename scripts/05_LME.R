library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)


# ============================================= #
# Data preparation                           ####
# ============================================= #

# Select and preprocess data from df_sum_analyte
df_lme <- df_sum_analyte %>%
  select(c(sample_date, treatment, analyte_name, mean)) %>%
  ungroup()

# Set treatment as factor with explicit levels
df_lme <- df_lme %>%
  mutate(treatment = factor(treatment, levels = c("WWTP-E", "O3", "AO", "CMF", "GAC")))

# Filter for measured values
df_lme <- df_lme %>% filter(mean > 0)

# ============================================= #
# Build Model                                ####
# ============================================= #

# Base model: mean ~ treatment + (1 | analyte_name)
lmer_chems <- lmer(mean ~ treatment + (1 | analyte_name), data = df_lme)
isSingular(lmer_chems)  # Check for singular fit
summary(lmer_chems)     # Model summary

# Extended model: add sample_date as fixed effect
lmer_chems_with_date <- lmer(
  mean ~ treatment + sample_date + (1 | analyte_name),
  data = df_lme
)
isSingular(lmer_chems_with_date)  

# Compare models
anova(lmer_chems, lmer_chems_with_date) # p > 0.05; sample_date does not improve model fit

# Summaries for both models
summary(lmer_chems)
summary(lmer_chems_with_date)

# ============================================= #
# Model diagnostics                          ####
# ============================================= #

# Leverage and Cook's distance plots
hat_values <- hatvalues(lmer_chems)
plot(hat_values, main = "Leverage Plot")
cooksd <- cooks.distance(lmer_chems)
plot(cooksd, main = "Cook's Distance Plot")

# Remove high-influence points (Cook's distance > 4/n)
df_reduced <- df_lme[-which(cooksd > 4/nrow(df_lme)), ]
lmer_reduced <- lmer(mean ~ treatment + (1 | analyte_name), data = df_reduced)
summary(lmer_reduced)

# Remove specific high-leverage analytes
df_without_high_leverage <- df_lme %>% filter(!analyte_name %in% c("Hydrochlorothiazide", "Benzotriazole", "Metoprolol"))
lmer_no_high_leverage <- lmer(mean ~ treatment + (1 | analyte_name), data = df_without_high_leverage)
summary(lmer_no_high_leverage)

# ============================================= #
# Random effect analysis                     ####
# ============================================= #

ranef(lmer_chems)   # Random effects estimates
VarCorr(lmer_chems) # Variance components

# Extract variance components as a data frame
var_components <- as.data.frame(VarCorr(lmer_chems)) %>% select(-c(var1, var2))
analyte_variance <- var_components$vcov[var_components$grp == "analyte_name"]
residual_variance <- var_components$vcov[var_components$grp == "Residual"]

# Analyte-to-analyte variability is large (~42% of residual variance)
explained_var <- round(analyte_variance / (analyte_variance + residual_variance) * 100)
unexplained_var <- 100 - explained_var

# ============================================= #
# Fixed effect analysis                      ####
# ============================================= #

anova(lmer_chems) # p < 0.05; at least one treatment differs significantly

# Pairwise comparisons
pair_comp_emmeans <- emmeans(lmer_chems, pairwise ~ treatment, adjust = "tukey")[[1]]  # Average concentrations
pair_comp_contrasts <- emmeans(lmer_chems, pairwise ~ treatment, adjust = "tukey")[[2]]  # Pairwise p-values
# Result: WWTP-E differs from all others (p < 0.05);
# no difference among O3, AO, CMF, GAC

# Pairwise comparisons without high leverage analytes
pair_comp_contrasts_no_lev <- emmeans(lmer_no_high_leverage, pairwise ~ treatment, adjust = "tukey")[[2]] 

# ============================================= #
# Model assumptions                          ####
# ============================================= #

# Normality and homoscedasticity checks
qqnorm(residuals(lmer_chems))
qqline(residuals(lmer_chems))
plot(fitted(lmer_chems), residuals(lmer_chems))
abline(h = 0)

# Log-transformation attempt
mod_log <- lmer(
  log10(mean) ~ treatment + (1 | analyte_name),
  data = df_lme
)
qqnorm(residuals(mod_log))
qqline(residuals(mod_log))
plot(fitted(mod_log), residuals(mod_log))
abline(h = 0)


# Compare log-transformed vs. original model
# Result: Same qualitative outcome; log model fits better (lower AIC/BIC);
# original scale retained for interpretability
anova(mod_log, lmer_chems) 
AIC(lmer_chems, mod_log) 
BIC(lmer_chems, mod_log) 
emmeans(mod_log, pairwise ~ treatment)  # Pairwise comparisons on log scale

