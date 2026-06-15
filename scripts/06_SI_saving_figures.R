library(patchwork)

# Note: Functions for saving images (i.e. save_image()) or 
# data frames (i.e save_df_csv) have been disabled by adding
# a comment symbol (#) so that # users can decide for themselves
# whether to save these objects.

# ================================================ #
# Figures and Tables: Main Manuscript           ####
# ================================================ #

## Section 3.1 Target Analysis – Measured OMP Concentrations and Removal  ####

### Figure 1 ####
p1_SI + p2_SI +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.key.size = unit(1, "cm")
    )
  )
# save_image("figure_1.png", width = 13, height = 8.5)

### Table 2 ####
target_data <- df_rem_total %>% rename("date" = "sample_date")

table_2_pre <- df_feat_n_wide %>% mutate(treatment = ifelse(treatment == "WWTP_E", "WWTP-E", as.character(treatment))) %>% 
  left_join(target_data, ., by = c("date", "treatment"))
table_2 <- table_2_pre %>% rename(
  Date = date,
  Treatment = treatment,
  `Summed concentration` = sum_conc,
  `Removal (%)` = rel_removal_total,
  Total = Full_NTS,
  Suspects = Suspect,
  TPs = TP
)
# save_df_csv(table_2, "table_2.csv")


## Section 3.2 Non-target analysis – Detected and Removed Features #### 

### Figure 2 ####
p1_trends + p2_trends +
  plot_layout(axis_titles = "collect",
              guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.key.size = unit(1.5, "cm")
    )
  ) 
# save_image("figure_2.png", width = 12.69, height = 8.27)


## Section 3.3 Bioassays – Detected Effects ####

### Figure 3 ####
wrap_plots(p1_BEQ + p2_BEQ_Eff) +
  plot_annotation(
    tag_levels = list(c("A", "B", "", "", "", "", "")),
    theme = theme(
      legend.position = "bottom",
      legend.box.spacing = unit(-0.5, "cm")
    )
  ) +
  plot_layout(guides = "collect")
# save_image("figure_3.png", width = 18, height = 9.5)

## Section 3.4 Potential Risks Based on Target Analysis ####

### Table 3 ####

tbl03_manuscript <- tbl03 %>% 
mutate(mean_sd = paste0(mean.RQ, " ± ", sd.RQ)) %>% 
  select(-matches("\\.RQ")) %>% 
  relocate(mean_sd, .after = sum_RQ) %>% 
  rename(Treatment = treatment,
         `Total RQ` = sum_RQ,
         `Mean ± SD` = mean_sd)
tbl03_manuscript 
# save_df_csv(tbl03_manuscript, "table_3.csv")


## Section 3.5 Potential Risks Based on Bioassays ####

## Figure 4 ####
p1_point_EBT_SW + p2_point_EBT_DW +
  plot_annotation(tag_levels = "A",
                  theme = theme(
                    legend.position = "bottom",
                    legend.box.spacing = unit(-0.1, "cm")
                    )
  ) +
  plot_layout(axis_titles = "collect", guides = "collect")
# save_image("Figure_4.png")

# =================================================== #
# Figures and Tables: Supplementary Information    ####
# =================================================== #

## Table S2.2 ####
# save_df_csv(tabSI_SW, "table_S2.2.csv")

## Table S2.3 ####
# save_df_csv(tabSI_DW, "table_S2.3.csv")

## Table S2.6 ####

tbl_S2.6 <- detection_overall %>% select(-3) %>% 
  mutate(detected_days = case_when(
    fod_detection == 6/6 ~ "6 of 6 days",
    fod_detection == 5/6 ~ "5 of 6 days",
    fod_detection == 4/6 ~ "4 of 6 days",
    fod_detection == 3/6 ~ "3 of 6 days",
    fod_detection == 2/6 ~ "2 of 6 days",
    fod_detection == 1/6 ~ "1 of 6 days"
  ))
# save_df_csv(tbl_S2.6, "table_S2.6.csv")

## Table S2.7 ####
tbl_S2.7 <- detection_treatment %>% select(-4) %>% 
  mutate(detected_days = case_when(
    
    # For WWTP-E, CMF
    fod_detection == 6/6 & treatment %in% c("WWTP-E", "CMF") ~ "6 of 6 days",
    fod_detection == 5/6 & treatment %in% c("WWTP-E", "CMF") ~ "5 of 6 days",
    fod_detection == 4/6 & treatment %in% c("WWTP-E", "CMF") ~ "4 of 6 days",
    fod_detection == 3/6 & treatment %in% c("WWTP-E", "CMF") ~ "3 of 6 days",
    fod_detection == 2/6 & treatment %in% c("WWTP-E", "CMF") ~ "2 of 6 days",
    fod_detection == 1/6 & treatment %in% c("WWTP-E", "CMF") ~ "1 of 6 days",
    
    # For O3
    fod_detection == 4/4 & treatment == "O3" ~ "4 of 4 days",
    fod_detection == 3/4 & treatment == "O3" ~ "3 of 4 days",
    fod_detection == 2/4 & treatment == "O3" ~ "2 of 4 days",
    fod_detection == 1/4 & treatment == "O3" ~ "1 of 4 days",
    fod_detection == 0/4 & treatment == "O3" ~ "0 of 4 days",
    
    # For AO
    fod_detection == 2/2 & treatment == "AO" ~ "2 of 2 days",
    fod_detection == 1/2 & treatment == "AO" ~ "1 of 2 days",
    fod_detection == 0/2 & treatment == "AO" ~ "0 of 2 days",
    
    # For GAC
    fod_detection == 5/5 & treatment == "GAC" ~ "5 of 6 days",
    fod_detection == 4/5 & treatment == "GAC" ~ "4 of 6 days",
    fod_detection == 3/5 & treatment == "GAC" ~ "3 of 6 days",
    fod_detection == 2/5 & treatment == "GAC" ~ "2 of 6 days",
    fod_detection == 1/5 & treatment == "GAC" ~ "1 of 6 days",
    fod_detection == 0/5 & treatment == "GAC" ~ "0 of 6 days"
  ))
tbl_S2.7
# save_df_csv(tbl_S2.7, "table_S2.7.csv")


## Table S2.8
tbl_S2.8 <- df_nts %>% filter(group %in% df_fgroups_rem_O3AO) %>% 
  distinct(group, .keep_all = TRUE) %>% 
  select(group, susp_name) %>% 
  left_join(., df_general_trend_long, by = "group") %>% 
  distinct(group, .keep_all = TRUE) %>% 
  select(group, ret, mz, susp_name, set)
tbl_S2.8
# save_df_csv(tbl_S2.8, "table_S2.8.csv")

## 3.1 Non-target data processing ####

### Figure S2 ####
p1_total_features_final + p2_avg_features_final + p3_data_size_final +
  plot_layout(guides = "collect", axis_titles = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"),
                  tag_levels = "A")
# save_image("figure_s2.png")

## 3.2 Target Analysis – Observed Removal  ####

### Table S11 ####
# Note: lmer_chems was generated in script 05
tbl_S11_prep <- summary(lmer_chems)$coefficients %>% as.data.frame()
tbl_S11 <- tbl_s11_prep %>% rownames_to_column(var = "Treatment") %>%
  rename(`p-value` = `Pr(>|t|)`,
         `t-value` = "t value") %>%
  mutate(
    `p-value` = ifelse(`p-value` < 0.001, "<0.001", round(`p-value`, 2)),
    across(c(Estimate, `Std. Error`, df, `t-value`), ~ round(., 2))
  )
tbl_S11
# save_df_csv(tbl_S11, "table_s11.csv")

  
### Table S12 ####
# Note: var_components, explained_var, and unexplained_var were generated in script 05
tbl_S12 <- var_components %>%
  rename(Group = grp,
         Variance = vcov,
         `Std. Dev.` = sdcor) %>%
  mutate(`Total Variance [%]` = c(explained_var, unexplained_var),
         across(c(Variance, `Std. Dev.`), ~ round(., 0)))
tbl_S12
# save_df_csv(tbl_S12, "table_s12.csv")


### Table S13 ####
# Note: pair_comp_contrasts and pair_comp_contrasts_no_lev were generated in script 05
tbl_S13 <- pair_comp_contrasts %>% as.data.frame() %>%
  rename(
    Contrast = contrast,
    Estimate = estimate,
    df = df,
    `t-value` = t.ratio,
    `p-value` = p.value
  ) %>%
  mutate(
    across(c(Estimate, SE, df, `t-value`), ~ round(., 2)),
    `p-value` = round(as.numeric(`p-value`), digits = 3),
    `p-value` = ifelse(`p-value` == 0.000, "<0.001", `p-value`)
    )
tbl_S13
# save_df_csv(tbl_S13, "table_s13.csv")

### Table S14 ####
tbl_S14 <- pair_comp_contrasts_no_lev %>% as.data.frame() %>%
  rename(
    Contrast = contrast,
    Estimate = estimate,
    df = df,
    `t-value` = t.ratio,
    `p-value` = p.value
  ) %>%
  mutate(
    across(c(Estimate, SE, df, `t-value`), ~ round(., 2)),
    `p-value` = round(as.numeric(`p-value`), digits = 3),
    `p-value` = ifelse(`p-value` == 0.000, "<0.001", `p-value`)
  )
tbl_S14
# save_df_csv(tbl_S14, "table_s14.csv")
  

### Figure S3 ####
figure_S3 + theme_bw(base_size = 22) +  
  theme(
    axis.text.y = element_text(size = 16))
# save_image("figure_S3.png", width = 8, height = 10)



## 3.3	 Non-target analysis – Detected and Removed Features ####

### Figure S4 ####
fig_S4 <- p_feat_n + p_feat_i +
  plot_annotation(tag_levels = "A",
                  theme = theme(legend.position = "bottom")) +
  plot_layout(guides = "collect")
fig_S4
# save_image("figure_S4.png", orientation = "landscape")

### Figure S5 ####
p_NTS_histo
# save_image("figure_S5.png")

### Figure S6 ####
p3_trends
# save_image("figure_S6.png")

### Figure S7 ####
p_NTS_boxplot
# save_image("figure_S7.png")

## 3.4	 Potential Risks Based on Bioassays ####

### Figure S8 #####
pDF10SW + pDF10DW + plot_annotation(tag_levels = "A",
                                    theme = theme(legend.position = "bottom")
) +
  plot_layout(axis_titles = "collect", guides = "collect")
# save_image("figure_S8.png")

## 3.5	 Implications for Water Reuse ####

### Figure S9 ####
OMP_pesticides
# save_image("figure_S9.png")
