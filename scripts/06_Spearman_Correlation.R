library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(forcats)
library(readxl)
library(rstatix)

# ## ============================================= #
# # Main Manuscript                             ####
# ## ============================================= #

# df_EBM_mean_eff was generated in script 03_EBM
# df_mean_intensity_spearman was generated in script 02_NTS_Analysis

df_EBM_mean_eff_spearman <- df_EBM_mean_eff %>% mutate(treatment = as.character(treatment))
df_mean_intensity_spearman <- df_mean_intensity_spearman %>% mutate(treatment = if_else(treatment == "WWTP_E", "WWTP-E", treatment))

df_CALUX_PXR <- df_EBM_mean_eff_spearman %>%
  filter(calux_assay == "PXR") %>%
  select(treatment, calux_assay, date, mean_eff) %>%
  ungroup()
df_CALUX_Nrf2 <- df_EBM_mean_eff_spearman %>%
  filter(calux_assay == "Nrf2") %>%
  select(treatment, calux_assay, date, mean_eff) %>%
  ungroup()

df_CALUX_NTS_PXR <- inner_join(df_mean_intensity_spearman, df_CALUX_PXR, by = c("date", "treatment"))
df_CALUX_NTS_Nrf <- inner_join(df_mean_intensity_spearman, df_CALUX_Nrf2, by = c("date", "treatment"))
df_CALUX_NTS <- rbind(df_CALUX_NTS_PXR, df_CALUX_NTS_Nrf)

treatments <- c("WWTP-E", "O3", "AO", "CMF", "GAC")
df_CALUX_NTS <- df_CALUX_NTS %>%
  mutate(treatment = factor(
    treatment,
    levels = treatments
  )) %>%
  select(-assigned_pol_mean)


n_pairs_all <- df_CALUX_NTS %>%
  group_by(group, calux_assay) %>%
  summarise(
    # Number of paired observations used by the Spearman test
    n_test = sum(
      is.finite(intensity) &
        is.finite(mean_eff)
    ),

    # Number of pairs where intensity and effects are not 0
    n_pairs_nonzero = sum(
      intensity != 0 &
        mean_eff != 0
    ),

    # Number of pairs where intensity is 0
    n_pairs_feature_zero = sum(
      mean_eff != 0 &
        intensity == 0
    ),
    .groups = "drop"
  ) %>%
  mutate(prop_nonzero = n_pairs_nonzero / n_test)
n_pairs_all

spearman <- df_CALUX_NTS %>%
  group_by(group, calux_assay) %>%
  rstatix::cor_test(
    intensity, mean_eff,
    method = "spearman"
  ) %>%
  group_by(calux_assay) %>%
  mutate(
    p_adj_BH = p.adjust(p, method = "BH")
  ) %>%
  arrange(p_adj_BH, desc(abs(cor)))

spearman_pairs <- spearman %>%
  left_join(
    n_pairs_all,
    by = c("group", "calux_assay")
  ) %>%
  select(
    group,
    n_pairs_nonzero,
    prop_nonzero,
    cor,
    statistic,
    p,
    p_adj_BH
  )

spearman_pairs %>%
  filter(
    calux_assay == "Nrf2",
    abs(cor) >= 0.8, p_adj_BH < 0.05
  ) # No feature below adjusted p-value

fgroups_PXR <- spearman_pairs %>%
  filter(
    calux_assay == "PXR",
    abs(cor) >= 0.7, p_adj_BH < 0.05
  ) %>%
  select(c(group, n_pairs_nonzero, cor, p_adj_BH, prop_nonzero)) %>%
  arrange(desc(abs(cor)))
fgroups_PXR


fgroups_PXR_cor_pos <- fgroups_PXR %>% filter(cor > 0)

fgroups_PXR_cor_neg <- fgroups_PXR %>% filter(cor < 0)

# Positive Correlation - Potentially effect inducing
p_PXR_positive_cor <- df_CALUX_NTS %>%
  filter(group %in% fgroups_PXR_cor_pos$group) %>%
  ggplot(aes(x = intensity, y = mean_eff / 1000, col = treatment)) +
  geom_jitter(height = 0.2, width = 2, alpha = 0.7) +
  scale_color_brewer(palette = "Set2") +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6)) +
  labs(
    x = "Intensity (×10⁶)",
    y = expression(paste("BEQ"["PXR"], " (µgL"^-1, ")")),
    col = NULL
  ) +
  theme(
    legend.position = "bottom",
    legend.key = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
p_PXR_positive_cor


# Negative Correlation - Potentially not effect inducing or potential TPs of effect inducing features
p_PXR_negative_cor <- df_CALUX_NTS %>%
  filter(group %in% fgroups_PXR_cor_neg$group) %>%
  ggplot(aes(x = intensity, y = mean_eff / 1000, col = treatment)) +
  geom_jitter(height = 0.2, width = 2, alpha = 0.7) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6)) +
  scale_color_brewer(palette = "Set2") +
  labs(
    x = "Intensity (×10⁶)",
    y = expression(paste("BEQ"["PXR"], " (µgL"^-1, ")")),
    col = NULL
  ) +
  theme(
    legend.position = "bottom",
    legend.key = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
p_PXR_negative_cor

# Extract feature groups from df_general_trend_long, seperate by set
fgroups_PXR_cor_pos_tbl <- df_general_trend_long %>%
  filter(group %in% fgroups_PXR_cor_pos$group) %>%
  distinct(group, .keep_all = TRUE) %>%
  select(group, ret, mz, susp_name, set) %>%
  mutate(correlation = "positive (rho > 0.7)")
fgroups_PXR_cor_neg_tbl <- df_general_trend_long %>%
  filter(group %in% fgroups_PXR_cor_neg$group) %>%
  distinct(group, .keep_all = TRUE) %>%
  select(group, ret, mz, susp_name, set) %>%
  mutate(correlation = "negative (rho < -0.7)")

tbl_S2_11 <- rbind(fgroups_PXR_cor_pos_tbl, fgroups_PXR_cor_neg_tbl)


# Plots ####
top_positive <- fgroups_PXR %>%
  slice_max(
    order_by = cor,
    n = 9,
    with_ties = FALSE
  ) %>%
  pull(group)

fig_S13 <- df_CALUX_NTS %>%
  filter(
    calux_assay == "PXR",
    group %in% top_positive
  ) %>%
  ggplot(aes(
    x = intensity,
    y = mean_eff / 1000,
    colour = treatment
  )) +
  geom_point(alpha = 0.7) +
  scale_color_brewer(palette = "Set2") +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6)) +
  facet_wrap(~group, scales = "free") +
  labs(
    x = "Intensity (×10⁶)",
    y = expression(paste("BEQ"["PXR"], " (µgL"^-1, ")")),
    colour = NULL
  ) +
  theme(
    legend.position = "bottom",
    legend.key = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1),
    strip.text.x = element_text(size = 14, face = "bold")
  )
fig_S13
