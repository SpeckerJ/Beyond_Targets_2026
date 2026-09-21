library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(readxl)
library(outliers)


## ============================================= #
## Read in data                               ####
## ============================================= #

# List all files in the  directory, read in data
raw_data_folder <- list.files("raw_data/EBM", full.names = TRUE)
raw_data_file_EBM <- raw_data_folder[str_detect(raw_data_folder, "raw_data_EBM.xlsx")]
raw_data_EBM <- read_excel(raw_data_file_EBM)

# assign to working df, convert to lowercase
df_raw_EBM <- raw_data_EBM %>% rename_with(str_to_lower)

## ============================================= #
## Data Wrangling                             ####
## ============================================= #

# Lookup tables for dates
lookup_dates <- data.frame(
  day = c(paste0("D", 1:6)),
  date = paste0(c("19.07.", "26.07.", "16.08.", "26.08.", "15.09.", "19.10."), "2022")
)

# Lookup table for EBT (Effect-Based Trigger) values
lookup_ebt <- df_raw_EBM %>%
  select(matches("calux|ebt")) %>%
  distinct() %>%
  mutate(
    ebt_sw_df10 = ebt_sw_ngl * 10,
    ebt_dw_df10 = ebt_dw_ngl * 10,
  )


# Convert results to numeric, flag LOQs, uniform unit
df_EBM <- df_raw_EBM %>%
  mutate(
    below_loq = ifelse(str_detect(result, "LOQ|<\\s*\\d"), 1, 0),
    result = as.numeric(result),
    result = ifelse(is.na(result), 0, result),
    result = ifelse(unit == "ug/l", result * 1000, result),
    loq = ifelse(unit == "ug/l", loq * 1000, loq),
    unit = ifelse(unit == "ug/l", "ng/l", unit),
    loq2 = loq/2,
    result_loq = as.numeric(ifelse(below_loq == 1, loq, result)),
    result_loq2 = as.numeric(ifelse(below_loq == 1, loq2, result))
  ) %>%
  relocate(loq2, .after = loq) %>% 
  relocate(result, .before = result_loq)
  

# Define order and change sample/treatment names
treatments <- c("WWTP-E", "O3", "AO", "CMF", "GAC") 
treat_sample_order <- c(treatments, "Field_B", "Lab_B_1L", "Lab_B_0.5L", "Pos_C")
df_EBM <- df_EBM %>%
  rename(treatment = sample) %>%
  mutate(
    treatment = factor(
      treatment,
      levels = treat_sample_order
    ),
    date = factor(date, levels = c(lookup_dates$date))
  )

# Generate df for blank samples
df_blanks_effects <- df_EBM %>%
  filter(str_detect(treatment, "_B")) %>%
  group_by(calux_assay) %>%
  summarise(
    mean_eff_blank = mean(result, na.rm = TRUE),
    mean_eff_blank_loq = mean(result_loq, na.rm = TRUE),
    mean_eff_blank_loq2 = mean(result_loq2, na.rm = TRUE)
  )

## ============================================= #
## Data Analysis                              ####
## ============================================= #

### Number of treatment samples < LOQ ####
df_EBM %>% filter(treatment %in% treatments) %>%
  filter(below_loq == 1) %>% 
  group_by(calux_assay, treatment, date) %>% 
  summarise(n_loq = sum(below_loq)) %>% 
  ungroup() %>% 
  mutate(n_loq_total = sum(n_loq))


### ============================================= #
### PFAS Outlier Analysis                      ####
### ============================================= #

# Generate df for PFAS-CALUX
df_PFAS <- df_EBM %>% filter(calux_assay == "PFAS")

# Data from the present study only and combined with de Schepper et al. (see manuscript)
df_WWTP_Specker <- df_PFAS %>% filter(treatment == "WWTP-E") %>% select(replicate, result)
de_schepper_2023 <- data.frame(replicate = c(NA,NA), result = c(8400, 4000))
df_WWTP_Specker_schepper <- rbind(df_WWTP_Specker, de_schepper_2023)

# Apply Dixon and Grubbs tests for outliers
dixon.test(df_WWTP_Specker$result, two.sided = TRUE) # p > 0.05; no outlier
grubbs.test(df_WWTP_Specker$result, two.sided = TRUE) # p > 0.05; no outlier
dixon.test(df_WWTP_Specker_schepper$result, two.sided = TRUE) # p < 0.05; outlier 
grubbs.test(df_WWTP_Specker_schepper$result, two.sided = TRUE) # p < 0.05; outlier 

# Remove outlier PFAS-CALUX
outlier_pos <- df_EBM %>%
  filter(calux_assay == "PFAS") %>%
  slice_max(result_loq, n = 1, with_ties = FALSE)
df_EBM <- df_EBM %>% filter(rowid != outlier_pos$rowid)


### ============================================= #
### Mean Effects                               ####
### ============================================= #

### ============================================= #
#### Effect per bioassay, day, treatment       ####
### ============================================= #

df_EBM_mean_eff <- df_EBM %>%
  filter(treatment %in% treatments) %>% 
  group_by(treatment, calux_assay, date) %>%
  summarise(
    mean_eff = mean(result, na.rm = TRUE),
    sd_eff = sd(result, na.rm = TRUE),
    
    mean_loq_eff = mean(result_loq, na.rm = TRUE),
    sd_loq_eff = sd(result_loq, na.rm = TRUE),
    
    mean_loq2_eff = mean(result_loq2, na.rm = TRUE),
    sd_loq2_eff = sd(result_loq2, na.rm = TRUE)
    
  ) %>%
  mutate(unit = "ng/L") %>%
  left_join(lookup_ebt, join_by(calux_assay))

### ============================================= #
### PFAS Effects with(out) outlier
### ============================================= #

df_removel_PFAS_with_outlier <- df_PFAS %>%
  filter(treatment %in% treatments) %>% 
  group_by(treatment, date) %>%
  summarise(
    mean_Eff = mean(result, na.rm = TRUE),
    sd_eff_PFAS_with_out = sd(result, na.rm = TRUE),
  ) %>%
  ungroup() %>% 
  mutate(unit = "ng/L",
         outlier = "outlier",
         baseline            = mean_Eff[treatment == "WWTP-E"],
         rel_removal_overall = (baseline - mean_Eff) / baseline * 100,
         rel_removal_overall = ifelse(rel_removal_overall == 0, NA, rel_removal_overall),
         rel_removal_overall = round(rel_removal_overall)) %>%
  arrange(date, treatment) %>% 
  pivot_longer(cols = c(rel_removal_overall, mean_Eff), names_to = "variable", values_to = "value") %>% 
  select(c(1,2, outlier, variable, value))

df_removal_PFAS_without_outlier <- df_PFAS %>%
  filter(rowid != outlier_pos$rowid, treatment %in% treatments) %>%
  group_by(treatment, date) %>%
  summarise(
    mean_Eff = mean(result, na.rm = TRUE),
    sd_eff_PFAS_no_out = sd(result, na.rm = TRUE),
  ) %>%
  ungroup() %>% 
  mutate(unit = "ng/L",
         outlier = "no_outlier",
         baseline            = mean_Eff[treatment == "WWTP-E"],
         rel_removal_overall = (baseline - mean_Eff) / baseline * 100,
         rel_removal_overall = ifelse(rel_removal_overall == 0, NA, rel_removal_overall),
         rel_removal_overall = round(rel_removal_overall)) %>%
  arrange(date, treatment) %>% pivot_longer(cols = c(rel_removal_overall, mean_Eff), names_to = "variable", values_to = "value") %>% 
  select(c(1,2, outlier, variable, value))
df_PFAS_rem <- rbind(df_removel_PFAS_with_outlier, df_removal_PFAS_without_outlier)

### ============================================= #
### Blank influence                            ####
### ============================================= #

df_EBM_mean_eff_blank <- left_join(df_EBM_mean_eff, df_blanks_effects, by = "calux_assay")

# Replace PAH BEQ after GAC by blank effect as BEQ Blank > BEQ sample
# (Effect is indistinguishable from blank; assume conservative worst-case)
# Note: SD values for PAH, GAC were not adjusted!
PAH_blank <- df_EBM_mean_eff_blank %>% filter(calux_assay == "PAH", treatment == "GAC") %>% pull(mean_eff_blank)
df_EBM_mean_eff_blank <- df_EBM_mean_eff_blank %>% mutate(
  mean_eff = case_when(
  calux_assay == "PAH" & treatment == "GAC" ~ PAH_blank,
  TRUE ~ mean_eff),
  mean_loq_eff = case_when(calux_assay == "PAH" & treatment == "GAC" ~ PAH_blank,
                         TRUE ~ mean_loq_eff),
  mean_loq2_eff = case_when(calux_assay == "PAH" & treatment == "GAC" ~ PAH_blank,
                            TRUE ~ mean_loq2_eff))

df_EBM_mean_eff_blank <- df_EBM_mean_eff_blank %>% mutate(
  mean_eff_blank_cor = mean_eff - mean_eff_blank,
  mean_loq_eff_blank_cor = mean_loq_eff - mean_eff_blank_loq,
  mean_loq2_eff_blank_cor = mean_loq2_eff - mean_eff_blank_loq2
  )

# Blank BEQ percentage of sample BEQ
df_blank_percentage <- df_EBM_mean_eff_blank %>% mutate(
  blank_per = round(mean_eff_blank/mean_eff * 100, digits = 2),
  blank_loq_per = round(mean_eff_blank_loq/mean_loq_eff * 100, digits = 2),
  blank_loq2_per = round(mean_eff_blank_loq2/mean_loq2_eff * 100, digits = 2),
) %>% select(c(1:3), matches("_per"))


### ============================================= #
### Effect Removal & EBT Assessment            ####
### ============================================= #

# Generate a list of all assays for further processing
df_list <- df_EBM_mean_eff_blank %>%
  filter(!is.na(treatment)) %>% # filter NAs due to FB, Lab Blanks etc.
  group_by(calux_assay) %>%
  group_split()


# Calculate mean removal per list element
df_list_output <- df_list %>%
  map(~ .x %>%
    group_by(date) %>%
    mutate(
      # Removal based on LOQ replacement with mean effects
      baseline            = mean_loq_eff[treatment == "WWTP-E"],
      rel_removal_overall = (baseline - mean_loq_eff) / baseline * 100,
      rel_removal_overall = ifelse(rel_removal_overall == 0, NA, rel_removal_overall),
      rel_removal_overall = round(rel_removal_overall),
                                  
      # Removal based on LOQ replacement with blank corrected effects
      baseline_blank            = mean_loq_eff_blank_cor[treatment == "WWTP-E"],
      rel_removal_overall_blank = (baseline_blank - mean_loq_eff_blank_cor) / baseline_blank * 100,
      rel_removal_overall_blank = ifelse(rel_removal_overall_blank == 0, NA, rel_removal_overall_blank),
      rel_removal_overall_blank = round(rel_removal_overall_blank)
    ) %>%
    arrange(date, treatment))

# Combine lists into one data frame
df_EBM_mean_eff <- bind_rows(df_list_output)

# Calculate mean, sd, and median removal
df_EBM_median_rem <- df_EBM_mean_eff %>%
  group_by(calux_assay, treatment, unit) %>%
  summarise(
    mean_rem = mean(rel_removal_overall, na.rm = TRUE),
    sd_rem = sd(rel_removal_overall, na.rm = TRUE),
    median_rem = median(rel_removal_overall, na.rm = TRUE),
    median_eff = median(mean_eff)
  )

# Calculate EBT exceedance
categorise_exceedance <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x <= 1   ~ "RC0",
    x <= 3   ~ "RC1",
    x <= 10  ~ "RC2",
    x <= 100 ~ "RC3",
    x > 100  ~ "RC4"
  )
}

df_EBM_mean_eff <- df_EBM_mean_eff %>%
  mutate(
    
    # Surface water exceedance
    ebt_sw_ex           = mean_eff      / ebt_sw_ngl,
    ebt_sw_ex_loq       = mean_loq_eff  / ebt_sw_ngl,
    ebt_sw_ex_loq2      = mean_loq2_eff / ebt_sw_ngl,
    ebt_sw_df10_ex      = mean_eff      / ebt_sw_df10,
    ebt_sw_df10_ex_loq  = mean_loq_eff  / ebt_sw_df10,
    ebt_sw_df10_ex_loq2 = mean_loq2_eff / ebt_sw_df10,
    
    # Drinking water exceedance
    ebt_dw_ex           = mean_eff      / ebt_dw_ngl,
    ebt_dw_ex_loq       = mean_loq_eff  / ebt_dw_ngl,
    ebt_dw_ex_loq2      = mean_loq2_eff / ebt_dw_ngl,
    ebt_dw_df10_ex      = mean_eff      / ebt_dw_df10,
    ebt_dw_df10_ex_loq  = mean_loq_eff  / ebt_dw_df10,
    ebt_dw_df10_ex_loq2 = mean_loq2_eff / ebt_dw_df10,
    
    # Surface water exceedance, blank corrected
    ebt_sw_ex_blank_cor           = mean_eff_blank_cor      / ebt_sw_ngl,
    ebt_sw_ex_loq_blank_cor       = mean_loq_eff_blank_cor  / ebt_sw_ngl,
    ebt_sw_ex_loq2_blank_cor      = mean_loq2_eff_blank_cor / ebt_sw_ngl,
    ebt_sw_df10_ex_blank_cor      = mean_eff_blank_cor      / ebt_sw_df10,
    ebt_sw_df10_ex_loq_blank_cor  = mean_loq_eff_blank_cor  / ebt_sw_df10,
    ebt_sw_df10_ex_loq2_blank_cor = mean_loq2_eff_blank_cor / ebt_sw_df10,
    
    # Drinking water exceedance, blank corrected
    ebt_dw_ex_blank_cor           = mean_eff_blank_cor      / ebt_dw_ngl,
    ebt_dw_ex_loq_blank_cor       = mean_loq_eff_blank_cor  / ebt_dw_ngl,
    ebt_dw_ex_loq2_blank_cor      = mean_loq2_eff_blank_cor / ebt_dw_ngl,
    ebt_dw_df10_ex_blank_cor      = mean_eff_blank_cor      / ebt_dw_df10,
    ebt_dw_df10_ex_loq_blank_cor  = mean_loq_eff_blank_cor  / ebt_dw_df10,
    ebt_dw_df10_ex_loq2_blank_cor = mean_loq2_eff_blank_cor / ebt_dw_df10,
    
    # Categorise all exceedance columns at once
    across(
      .cols = matches("^ebt_(sw|dw)(_df10)?_ex(_loq2?)?(_blank_cor)?$"),
      .fns   = categorise_exceedance,
      .names = "{.col}_cat"
    )
  )

### ============================================= #
### LOQ and blank influence                    ####
### ============================================= #

# Influence of replacing values <LOQ with the LOQ or LOQ/2
tbl_LOQ <- df_EBM_mean_eff %>%
  select(treatment, date, calux_assay, matches("_ex|ex_")) %>%
  select(-ends_with("_cat"),
         -contains("_blank_")) %>%
  pivot_longer(
    cols = contains("ex"),
    names_to = "variable",
    values_to = "value"
  ) %>%
  left_join(
    df_EBM_mean_eff %>%
      select(treatment, date, calux_assay, ends_with("_cat")) %>%
      pivot_longer(
        cols = ends_with("_cat"),
        names_to = "variable_cat",
        values_to = "category"
      ) %>%
      mutate(variable = str_remove(variable_cat, "_cat")),
    by = c("treatment", "date", "calux_assay", "variable")
  ) %>%
  filter(
    !is.na(value),
   calux_assay %in% c("ERα", "PXR", "Nrf2"),
   treatment %in% c("CMF", "GAC")
  ) %>%
  mutate(
    value = round(value, digits = 2),
    matrix_df = paste(
      ifelse(str_detect(variable, "_sw_"), "SW", "DW"),
      ifelse(str_detect(variable, "_df10_"), "1", "0")
    )
  ) %>%
  select(-variable_cat)
tbl_LOQ_comp <- tbl_LOQ %>%
  filter(
    calux_assay %in% c("ERα", "Nrf2", "PXR"),
    !is.na(value)
  ) %>%
  group_by(treatment, date, calux_assay, matrix_df) %>%
  summarise(
    mean_value = round(mean(value), digits = 2),
    sd = round(sd(value), digits = 2),
    rsd = round((sd / mean_value) * 100, digits = 2),
    values = paste(value, collapse = "; "),
    n.categories = n_distinct(category, na.rm = TRUE),
    
    category = paste(
      sort(unique(category[!is.na(category)])),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>% select(-n.categories) %>% 
  filter(sd > 0)
tbl_LOQ_comp

tbl_blank <- df_EBM_mean_eff %>%
  select(treatment, date, calux_assay, matches("_ex|ex_")) %>%
  select(-ends_with("_cat")) %>%
  pivot_longer(
    cols = contains("ex"),
    names_to = "variable",
    values_to = "value"
  ) %>%
  left_join(
    df_EBM_mean_eff %>%
      select(treatment, date, calux_assay, ends_with("_cat")) %>%
      pivot_longer(
        cols = ends_with("_cat"),
        names_to = "variable_cat",
        values_to = "category"
      ) %>%
      mutate(variable = str_remove(variable_cat, "_cat")),
    by = c("treatment", "date", "calux_assay", "variable")
  ) %>%
  mutate(
    value = round(value, digits = 2),
    matrix_df = paste(
      ifelse(str_detect(variable, "_sw_"), "SW", "DW"),
      ifelse(str_detect(variable, "_df10_"), "1", "0"),
      ifelse(str_detect(variable, "_blank_"), "B", "noB")
    )
  ) %>%
  select(-variable_cat)

tbl_blank_comp <- tbl_blank %>%
  filter(
    calux_assay != "Cytotox",
    !is.na(value)
  ) %>%
  group_by(treatment, date, calux_assay, matrix_df) %>%
  summarise(
    mean_value = round(mean(value), digits = 2),
    sd = round(sd(value), digits = 2),
    rsd = round((sd / mean_value) * 100, digits = 2),
    values = paste(value, collapse = "; "),
    n.categories = n_distinct(category, na.rm = TRUE),
    
    category = paste(
      sort(category[!is.na(category)]),
      collapse = ", "
    ),
    
    .groups = "drop"
  ) %>% 
  filter(n.categories > 1) %>% select(-n.categories) %>% 
  mutate(matrix_df = str_remove(matrix_df, " B"))
tbl_blank_comp


# ============================================== #
# Plots                                       ####
# ============================================== #



## ============================================= #
### BEQ (Biological Equivalent Concentration) ####
## ============================================= #

# Change level (plot) order
df_EBM_mean_eff <- df_EBM_mean_eff %>%
  mutate(calux_assay = factor(calux_assay, levels = c(
    "PXR", "Nrf2", "Cytotox", "PAH", "ERα", "PFAS"
  )))

# Define levels
treatment_levels <- c("WWTP-E", "O3", "AO", "CMF", "GAC")

# Define date colors for plotting
date_levels <- levels(df_EBM_mean_eff$date)
date_cols <- brewer.pal(max(3, length(date_levels)), "Set2")[seq_along(date_levels)]
names(date_cols) <- date_levels
pos <- position_dodge(width = 0.3)
variable_to_plot <- "mean_loq_eff"

# Plot every subplot individually

#### ============================================ #
#### PXR                                       ####
#### ============================================ #

# Add dummy rows to PXR for 15.10.2022 to get a complete legend for all plots
dummy_row <- df_EBM_mean_eff %>% filter(calux_assay == "PFAS", treatment == "GAC")
dummy_row <- dummy_row %>% mutate(calux_assay = "PXR",
                                  across(where(is.numeric), ~ NA)
)
df_EBM_plot <- rbind(dummy_row, df_EBM_mean_eff)
df_EBM_plot <- df_EBM_plot %>% mutate(
  treatment = factor(treatment, levels = treatment_levels)) %>%
  complete(treatment, calux_assay)


p_BEQ_PXR <- df_EBM_plot %>%
  filter(!is.na(treatment), calux_assay == "PXR") %>%
  mutate(variable_to_plot = .data[[variable_to_plot]] / 1000) %>% # transform to µg/L
  ggplot(aes(x = treatment, y = variable_to_plot, fill = date, group = date, col = date)) +
  geom_point(
    size = 5,
    pch = 21,
    position = pos,
    na.rm = TRUE,
    alpha = 0.8,
    col = "black"
  ) +
  scale_fill_manual(values = date_cols, drop = FALSE) +
  scale_y_log10(limits = c(0.6, NA), breaks = c(1, 10, 100)) +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  geom_hline(aes(yintercept = 3,   linetype = "EBT_SW"), col = "lightsalmon", key_glyph = "path") +
  geom_hline(aes(yintercept = 0.7, linetype = "EBT_DW"), col = "lightblue",   key_glyph = "path") +
  scale_linetype_manual(
    name = NULL,
    values = c("EBT_SW" = "dashed", "EBT_DW" = "solid"),  # both dashed in plot
    labels = c(
      "EBT_SW" = expression(EBT[SW]),
      "EBT_DW" = expression(EBT[DW])
    )
  ) +
  guides(
    fill     = guide_legend(order = 1),
    linetype = guide_legend(
      order = 2,
      nrow = 2,
      override.aes = list(
        linetype  = c(1, 1),        # both solid in legend
        linewidth = c(1, 1),        # same size
        colour    = c("lightsalmon", "lightblue")  # retain line colours
      )
    ),
    colour = "none"
  ) +
  labs(
    x = NULL,
    y = expression(paste("µgL"^-1,"")),
    fill = NULL
  ) +
  guides(
    fill     = guide_legend(order = 1),
    linetype = guide_legend(
      order = 2,
      nrow = 2,
      override.aes = list(
        linetype  = c(1, 1),        # both solid in legend
        linewidth = c(1, 1),        # same size
        colour    = c("lightsalmon", "lightblue")  # retain line colours
      )
    ),
    colour = "none"
  ) +
  theme(
    legend.position = "bottom",
    legend.key = element_rect(fill = "white"),
    strip.text.x = element_text(face = "bold"),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey98")
  )  +
  facet_grid(. ~ calux_assay)
p_BEQ_PXR

#### ============================================ #
#### PAH                                       ####
#### ============================================ #
p_BEQ_PAH <- df_EBM_plot %>%
  filter(!is.na(treatment), calux_assay == "PAH") %>%
  mutate(variable_to_plot = .data[[variable_to_plot]] / 1000) %>% # transform to µg/L
  ggplot(aes(x = treatment, y = variable_to_plot, fill = date, group = date, col = date)) +
  geom_point(
    size = 5,
    pch = 21,
    position = pos,
    na.rm = TRUE,
    alpha = 0.8,
    col = "black"
  ) +
  scale_color_manual(values = date_cols, drop = FALSE) +
  scale_fill_manual(values = date_cols, drop = FALSE) +
  scale_y_log10(limits = c(5/1000, NA)) +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  geom_hline(yintercept = 6.2 / 1000, linetype = "dashed", col = "lightsalmon") + 
  geom_hline(yintercept = 19 / 1000, linetype = "solid", col = "lightblue") +
  labs(
    x = NULL,
    y = expression(paste("µgL"^-1, "")),
    fill = NULL
  ) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(face = "bold"),
    axis.text.x = element_text(angle = 40, hjust = 1),
    axis.title.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey98")
  )  +
  facet_grid(. ~ calux_assay)
p_BEQ_PAH

#### ============================================ #
#### Nrf2                                      ####
#### ============================================ #
p_BEQ_Nrf <- df_EBM_plot %>%
  filter(!is.na(treatment), calux_assay == "Nrf2") %>%
  mutate(variable_to_plot = .data[[variable_to_plot]] / 1000) %>% # transform to µg/L
  ggplot(aes(x = treatment, y = variable_to_plot, fill = date, group = date, col = date)) +
  geom_point(
    size = 5,
    pch = 21,
    position = pos,
    na.rm = TRUE,
    alpha = 0.8,
    col = "black"
  ) +
  scale_color_manual(values = date_cols, drop = FALSE) +
  scale_fill_manual(values = date_cols, drop = FALSE) +
  scale_y_log10(limits = c(9, NA)) +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  geom_hline(yintercept = 10, linetype = "dashed", col = "lightsalmon") +
  geom_hline(yintercept = 137, linetype = "solid", col = "lightblue") +
  labs(
    x = NULL,
    y = expression(paste("µgL"^-1, "")),
    fill = NULL
  ) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(face = "bold"),
    axis.text.x = element_text(angle = 40, hjust = 1),
    axis.title.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey98")
  )  +
  facet_grid(. ~ calux_assay)
p_BEQ_Nrf

#### ============================================ #
#### ERa                                       ####
#### ============================================ #
p_BEQ_ERa <- df_EBM_plot %>%
  filter(!is.na(treatment), calux_assay == "ERα") %>%
  mutate(variable_to_plot = .data[[variable_to_plot]]) %>% # Keep ng/L
  ggplot(aes(x = treatment, y = variable_to_plot, fill = date, group = date, col = date)) +
  geom_point(
    size = 5,
    pch = 21,
    position = pos,
    na.rm = TRUE,
    alpha = 0.8,
    col = "black"
  ) +
  scale_color_manual(values = date_cols, drop = FALSE) +
  scale_fill_manual(values = date_cols, drop = FALSE) +
  scale_y_log10(limits = c(0.0073, NA)) +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  geom_hline(yintercept = 0.1, linetype = "dashed", col = "lightsalmon") +
  geom_hline(yintercept = 0.0083, linetype = "solid", col = "lightblue") +
  labs(
    x = NULL,
    y = expression(paste("ngL"^-1, "")),
    fill = NULL
  ) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(face = "bold"),
    axis.text.x = element_text(angle = 40, hjust = 1),
    axis.title.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey98")
  )  +
  facet_grid(. ~ calux_assay)
p_BEQ_ERa

#### ============================================ #
#### Cytotox                                   ####
#### ============================================ #
p_BEQ_Cyt <- df_EBM_plot %>%
  filter(!is.na(treatment), calux_assay == "Cytotox") %>%
  mutate(variable_to_plot = .data[[variable_to_plot]]) %>% # transform to µg/L
  ggplot(aes(x = treatment, y = variable_to_plot, fill = date, group = date, col = date)) +
  geom_point(
    size = 5,
    pch = 21,
    position = pos,
    na.rm = TRUE,
    alpha = 0.8,
    col = "black"
  ) +
  # scale_color_manual(values = date_cols, drop = FALSE) +
  scale_fill_manual(values = date_cols, drop = FALSE) +
  scale_y_log10() +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  labs(
    x = NULL,
    y = expression(paste("µgL"^-1, "")),
    fill = NULL
  ) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(face = "bold"),
    axis.text.x = element_text(angle = 40, hjust = 1),
    axis.title.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey98")
  )  +
  facet_grid(. ~ calux_assay)
p_BEQ_Cyt

#### ============================================ #
#### PFAS                                      ####
#### ============================================ #
p_BEQ_PFAS <- df_EBM_plot %>%
  filter(!is.na(treatment), calux_assay == "PFAS") %>%
  mutate(variable_to_plot = .data[[variable_to_plot]] / 1000) %>% # transform to µg/L
  ggplot(aes(x = treatment, y = variable_to_plot, fill = date, group = date, col = date)) +
  geom_point(
    size = 5,
    pch = 21,
    position = pos,
    na.rm = TRUE,
    alpha = 0.8,
    col = "black"
  ) +
  scale_color_manual(values = date_cols, drop = FALSE) +
  scale_fill_manual(values = date_cols, drop = FALSE) +
  scale_y_log10(limits = c(0.61, NA)) +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  geom_hline(yintercept = 0.71, linetype = "dashed", col = "lightsalmon") +
  geom_hline(yintercept = 3, linetype = "solid", col = "lightblue") +      
  labs(
    x = NULL,
    y = expression(paste("µgL"^-1, "")),
    fill = NULL
  ) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(face = "bold"),
    axis.text.x = element_text(angle = 40, hjust = 1),
    axis.title.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey98")
  )  +
  facet_grid(. ~ calux_assay)
p_BEQ_PFAS

#### ============================================ #
#### Combine BEQ plots                         ####
#### ============================================ #

# Remove x-axis tick labels from PXR, Nrf2, and Cytotox-plots
p_BEQ_PXR_top <- p_BEQ_PXR +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )
p_BEQ_Nrf_top <- p_BEQ_Nrf +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )
p_BEQ_Cyt_top <- p_BEQ_Cyt +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

# Three rows, # = empty space
design_fig3_BEQ <- "
A#B#C
#####
D#E#F
"

subpannels_fig3 <- wrap_plots(
  p_BEQ_PXR_top,
  p_BEQ_Nrf_top,
  p_BEQ_Cyt_top,
  p_BEQ_PAH,
  p_BEQ_ERa,
  p_BEQ_PFAS,
  design = design_fig3_BEQ
) +
  plot_layout(
    heights = c(1, -0.06, 1),
    widths = c(1, -0.06, 1, -0.06, 1),
    guides = "collect",
    axis_titles = "collect_x"
  ) +
  plot_annotation(
    theme = theme(
      legend.position = "bottom",
      legend.box.spacing = unit(0, "pt"),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15),
      legend.key.size = unit(1, "cm"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "cm"),
      legend.margin = margin(0, 0, 0, 0, unit = "cm")
    )
  )

subpannels_fig3

blanklabel_BEQ <- ggplot() +
  labs(y = "BEQ") +
  theme_classic() +
  guides(x = "none", y = "none") +
  theme(axis.title = element_text(size = 26))
blanklabel_BEQ
fig_3A <- blanklabel_BEQ + subpannels_fig3 + plot_layout(widths=c(1,1000))
fig_3A

## ============================================= #
### Effect removal                            ####
## ============================================= #

# Generate every plot individually to control axis

#### ============================================ #
#### PXR                                       ####
#### ============================================ #
p_Eff_rem_PXR <- df_EBM_plot %>%
  filter(!is.na(treatment), treatment != "WWTP-E", calux_assay == "PXR") %>%
  mutate(variable_to_plot = .data[[variable_to_plot]] / 1000) %>% # transform to µg/L
  ggplot(aes(x = treatment, y = rel_removal_overall, group = date)) +
  geom_point(
    aes(fill = date),
    position = pos,
    size = 5,
    shape = 21,
    colour = "black",
    show.legend = TRUE
  ) +
  scale_fill_manual(values = date_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(0,100),
    breaks = c(seq(0, 100, by = 25))
  ) +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  labs(x = "", y = "Effect removal (%)", fill = "") +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    panel.background = element_rect(fill = "grey98")
  ) +
  facet_wrap(~calux_assay)
p_Eff_rem_PXR

#### ============================================ #
#### Nrf2                                      ####
#### ============================================ #
p_Eff_rem_Nrf <- df_EBM_plot %>%
  filter(!is.na(treatment), treatment != "WWTP-E", calux_assay == "Nrf2") %>%
  mutate(variable_to_plot = .data[[variable_to_plot]] / 1000) %>% # transform to µg/L
  ggplot(aes(x = treatment, y = rel_removal_overall, group = date)) +
  geom_point(
    aes(fill = date),
    position = pos,
    size = 5,
    shape = 21,
    colour = "black",
    show.legend = TRUE
  ) +
  scale_fill_manual(values = date_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(-150,100),
    breaks = c(seq(-150, 100, by = 50))
  ) +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  labs(x = "", y = "Effect removal (%)", fill = "") +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    panel.background = element_rect(fill = "grey98")
  ) +
  facet_wrap(. ~ calux_assay)
p_Eff_rem_Nrf

#### ============================================ #
#### Cytotox                                   ####
#### ============================================ #
p_Eff_rem_Cyto <- df_EBM_plot %>%
  filter(!is.na(treatment), treatment != "WWTP-E", calux_assay == "Cytotox") %>%
  mutate(variable_to_plot = .data[[variable_to_plot]] / 1000) %>% # transform to µg/L
  ggplot(aes(x = treatment, y = rel_removal_overall, group = date)) +
  geom_point(
    aes(fill = date),
    position = pos,
    size = 5,
    shape = 21,
    colour = "black",
    show.legend = TRUE
  ) +
  scale_fill_manual(values = date_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(-160,100),
    breaks = c(seq(-150, 100, by = 50))
  ) +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  labs(x = "", y = "Effect removal (%)", fill = "") +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    panel.background = element_rect(fill = "grey98")
  ) +
  facet_wrap( ~ calux_assay)
p_Eff_rem_Cyto

#### ============================================ #
#### PAH                                       ####
#### ============================================ #
p_Eff_rem_PAH <- df_EBM_plot %>%
  filter(!is.na(treatment), treatment != "WWTP-E", calux_assay == "PAH") %>%
  ggplot(aes(x = treatment, y = rel_removal_overall, group = date)) +
  geom_point(
    aes(fill = date),
    position = pos,
    size = 5,
    shape = 21,
    colour = "black",
    show.legend = TRUE
  ) +
  scale_fill_manual(values = date_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(70,100),
    breaks = c(seq(70, 100, by = 15))
  ) +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  labs(x = "", y = "Effect removal (%)", fill = "") +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    panel.background = element_rect(fill = "grey98")
  ) +
  facet_wrap(. ~ calux_assay)
p_Eff_rem_PAH

#### ============================================ #
#### ERa                                       ####
#### ============================================ #
p_Eff_rem_ERa <- df_EBM_plot %>%
  filter(!is.na(treatment), treatment != "WWTP-E", calux_assay == "ERα") %>%
  ggplot(aes(x = treatment, y = rel_removal_overall, group = date)) +
  geom_point(
    aes(fill = date),
    position = pos,
    size = 5,
    shape = 21,
    colour = "black",
    show.legend = TRUE
  ) +
  scale_fill_manual(values = date_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(70,100),
    breaks = c(seq(70, 100, by = 15))
  ) +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  labs(x = "", y = "Effect removal (%)", fill = "") +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    panel.background = element_rect(fill = "grey98")
  ) +
  facet_wrap(. ~ calux_assay)
p_Eff_rem_ERa

#### ============================================ #
#### PFAS                                      ####
#### ============================================ #
p_Eff_rem_PFAS <- df_EBM_plot %>%
  filter(!is.na(treatment), treatment != "WWTP-E", calux_assay == "PFAS") %>%
  ggplot(aes(x = treatment, y = rel_removal_overall, group = date)) +
  geom_point(
    aes(fill = date),
    position = pos,
    size = 5,
    shape = 21,
    colour = "black",
    show.legend = TRUE
  ) +
  scale_fill_manual(values = date_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(0,50),
    breaks = c(seq(0, 50, by = 25))
  ) +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  labs(x = "", y = "Effect removal (%)", fill = "") +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    panel.background = element_rect(fill = "grey98")
  ) +
  facet_wrap(. ~ calux_assay)
p_Eff_rem_PFAS

#### ============================================ #
#### Combine effect plots                      ####
#### ============================================ #

# Remove x-axis tick labels from PXR, Nrf2, and Cytotox-plots
p_Eff_rem_PXR_top <- p_Eff_rem_PXR +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )
p_Eff_rem_Nrf_top <- p_Eff_rem_Nrf +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )
p_Eff_rem_Cyto_top <- p_Eff_rem_Cyto +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

# Three rows, # = empty space
design_fig3_B <- "
A#B#C
#####
D#E#F
"
fig3_B <- wrap_plots(
  p_Eff_rem_PXR_top,
  p_Eff_rem_Nrf_top,
  p_Eff_rem_Cyto_top,
  p_Eff_rem_PAH,
  p_Eff_rem_ERa,
  p_Eff_rem_PFAS,
  design = design_fig3_B
) +
  plot_layout(
    widths = c(1, 0.03, 1, 0.03, 1),
    heights = c(1, 0.06, 1),
    guides = "collect",
    axis_titles = "collect_y"
  ) +
  plot_annotation(
    theme = theme(
      legend.position = "none",
      legend.box.spacing = unit(0, "pt"),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15),
      legend.key.size = unit(1, "cm"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "cm"),
      legend.margin = margin(0, 0, 0, 0, unit = "cm")
    )
  )
fig3_B



### ============================================= #
### PFAS - Outlier Comparison                  ####
### ============================================= #
my_cols <- brewer.pal(3, "Set1")[c(1, 2)] # red, blue


p1_PFAS_A <- df_PFAS_rem %>% 
  mutate(outlier_status = factor(
    outlier,
    levels = c("no_outlier", "outlier"),
    labels = c("Excluded", "Included"))
  ) %>% 
  filter(variable == "mean_Eff") %>% 
  ggplot(
    aes(
      x = treatment,
      y = value / 1000,
      colour = outlier_status,
      shape = outlier_status
    )
  ) +
  geom_point(size = 6,
             position = position_dodge(width = 0.2)
  ) +
  scale_colour_manual(
    values = c(
      "Excluded" = my_cols[2],
      "Included" = my_cols[1]
    ),
    name = "Outlier:"
  ) +
  scale_shape_manual(
    values = c(
      "Excluded" = 17,
      "Included" = 16
    ),
    name = "Outlier:"
  ) +
  labs(
    x = NULL,
    y = expression(paste("BEQ"["PFAS"], " (" * mu * "g L"^-1, ")"))
  ) +
  theme(
    text = element_text(size = 30),
    legend.key = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey98")
  )
p1_PFAS_A

p1_PFAS_B <- df_PFAS_rem %>% 
  mutate(outlier_status = factor(
    outlier,
    levels = c("no_outlier", "outlier"),
    labels = c("Excluded", "Included"))
  ) %>% 
  filter(variable == "rel_removal_overall", !is.na(value)) %>% 
  ggplot(
    aes(
      x = treatment,
      y = value,
      colour = outlier_status,
      shape = outlier_status
    )
  ) +
  geom_point(size = 6) +
  scale_colour_manual(
    values = c(
      "Excluded" = my_cols[2],
      "Included" = my_cols[1]
    ),
    name = "Observation"
  ) +
  scale_shape_manual(
    values = c(
      "Excluded" = 17,
      "Included" = 16
    ),
    name = "Observation"
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 25)
  ) +
  coord_cartesian(
    ylim = c(0, 100)
  ) +
  labs(
    x = NULL,
    y = "Effect removal (%)"
  ) +
  theme(
    text = element_text(size = 30),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey98")
  )
p1_PFAS_B



## ============================================= #
## Risk Analysis                              ####
## ============================================= #

## ============================================= #
## EBTs (Effect-Based Triggers)               ####
## ============================================= #


### ============================================= #
### Without dilution factor (DF)               ####
### ============================================= #

# Define dodging position
pos2 <- position_dodge2(width = 0.6)

# Define treatment levels and colors
levels_treats <- df_EBM_mean_eff %>% filter(treatment %in% c("WWTP-E", "AO", "O3", "CMF", "GAC"))
treat_levels <- factor(unique(levels_treats$treatment))
treat_cols <- brewer.pal(max(3, length(treat_levels)),
                             "Spectral")[seq_along(treat_levels)]
names(treat_cols) <- treat_levels


# Calculate midpoint for band labels
band_labels <- tibble::tribble(
  ~label, ~ymin, ~ymax,
  "C[0]", 0.45, 1,
  "C[1]", 1, 3,
  "C[2]", 3, 10,
  "C[3]", 10, 100,
  "C[4]", 100, 1000
)
band_labels <- band_labels %>%
  mutate(y = sqrt(ymin * ymax))


# Plot EBT exceedance for surface water
p1_point_EBT_SW <- df_EBM_mean_eff %>% 
  filter(!is.na(ebt_sw_ngl)) %>%
  mutate(calux_assay = factor(calux_assay, levels = c("PXR", "Nrf2", "PAH", "ERα", "PFAS"))) %>% 
  ggplot(aes(x = calux_assay, y = ebt_sw_ex)) +
  geom_hline(yintercept = c(1, 3, 10, 100), col = "grey50", linetype = "dashed") +
  geom_point(
    aes(fill = treatment),
    position = pos2,
    size = 6,
    shape = 21,
    colour = "black",
    alpha = 0.8
  ) +
  geom_text(
    data = band_labels,
    aes(x = Inf, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 1.1,       
    fontface = "bold",
    size = 10,
    parse = TRUE
  ) +
  coord_cartesian(
    clip = "off" # drawing outside panel
  ) +
  labs(
    x = NULL,
    y = "EBT-Exceedance",
    fill = NULL
  ) +
  scale_fill_manual(
    values = treat_cols, drop = TRUE,
    labels = c(
      "O3" = expression(paste("O"[3]))
    )
  ) +
  scale_y_log10(
    breaks = c(1, 3, 10, 100),
    limits = c(0.30, 320)
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey98")
  )
p1_point_EBT_SW

# Plot EBT exceedance for drinking water
p2_point_EBT_DW <- df_EBM_mean_eff %>%
  filter(!is.na(ebt_dw_ngl)) %>%
  mutate(calux_assay = factor(calux_assay, levels = c("PXR", "Nrf2", "PAH", "ERα", "PFAS"))) %>% 
  ggplot(aes(x = calux_assay, y = ebt_dw_ex)) +
  geom_hline(yintercept = c(1, 3, 10, 100), col = "grey50", linetype = "dashed") +
  geom_point(
    aes(fill = treatment),
    position = pos2,
    size = 6,
    shape = 21,
    colour = "black",
    alpha = 0.8
  ) +
  geom_text(
    data = band_labels,
    aes(x = Inf, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 1.1,       
    fontface = "bold",
    size = 10,
    parse = TRUE
  ) +
  coord_cartesian(
    clip = "off" # drawing outside panel
  ) +
  labs(
    x = NULL,
    y = "EBT-Exceedance",
    fill = NULL
  ) +
  scale_fill_manual(
    values = treat_cols, drop = TRUE,
    labels = c(
      "O3" = expression(paste("O"[3]))
    )
  ) +
  scale_y_log10(
    breaks = c(1, 3, 10, 100),
    limits = c(0.30, 320)
  ) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey98")
  )
p2_point_EBT_DW


### ============================================= #
### With DF = 10                               ####
### ============================================= #

# Define dodging position
pos3 <- position_dodge2(width = 1)

# Calculate upper point for band labels
# (below RC lines)
band_labels_DF10 <- tibble::tibble(
  label = c("RC[0]", "RC[1]", "RC[2]", "RC[3]", "RC[4]"),
  line  = c(1, 3, 10, 100, 320)
) %>%
  mutate(
    y = line / 1.03
  )
band_labels_DF10[5,3] <- Inf

# Prepare data for EBT exceedance with DF10
df_EBT_SW <- df_EBM_mean_eff %>%
  filter(!is.na(ebt_sw_ngl)) %>%
  select(treatment, calux_assay, ebt_sw_ex, ebt_sw_df10_ex) %>% 
  pivot_longer(cols = c(ebt_sw_ex, ebt_sw_df10_ex), names_to = "variable", values_to = "values") %>% 
  mutate(calux_assay = case_when(
    str_detect(variable, "_sw_df10") ~ paste0(calux_assay, "_DF10"),
    TRUE ~ calux_assay
  ))

df_EBT_DW <- df_EBM_mean_eff %>%
  filter(!is.na(ebt_dw_ngl)) %>%
  select(treatment, calux_assay, ebt_dw_ex, ebt_dw_df10_ex) %>% 
  pivot_longer(cols = c(ebt_dw_ex, ebt_dw_df10_ex), names_to = "variable", values_to = "values") %>% 
  mutate(calux_assay = case_when(
    str_detect(variable, "_dw_df10") ~ paste0(calux_assay, "_DF10"),
    TRUE ~ calux_assay
  ))

df_EBT_SW_DW_DF10 <- rbind(df_EBT_DW, df_EBT_SW)


treatment_levels <- levels(df_EBM_mean_eff$treatment)[1:5]
assay_levels_DF <- c("PXR", "PXR_DF10",
                     "Nrf2", "Nrf2_DF10",
                     "PAH","PAH_DF10",
                     "ERα", "ERα_DF10",
                     "PFAS", "PFAS_DF10")
treatment_shapes <- setNames(
  c(21, 22, 23, 24, 25),
  treatment_levels
)

# Add dummy rows for SW and DW for the 19.10.2022 to get a complete legend for all plots
dummy_rows <- data.frame(date = c("19.10.2022", "19.10.2022"), treatment = c("GAC", "GAC"), calux_assay = c("PXR", "PXR"),
                         variable = c("ebt_dw_ex", "ebt_sw_ex"), values = c(-1,-1))

df_EBT_DF_plot <- rbind(dummy_rows, df_EBT_SW_DW_DF10)
df_EBT_DF_plot <- df_EBT_DF_plot %>% mutate(
  treatment = factor(treatment, levels = treatment_levels)) %>%
  complete(treatment, calux_assay)

# Plot EBT exceedance for surface water
p1_point_EBT_SW_date <- df_EBT_DF_plot %>%
  mutate(calux_assay = factor(calux_assay, levels = assay_levels_DF)) %>% 
  filter(!is.na(values), str_detect(variable, "_sw_")) %>%
  ggplot(aes(x = calux_assay, y = values)
  ) +
  geom_hline(
    yintercept = c(1, 3, 10, 100),
    colour = "grey50",
    linetype = "dashed"
  ) +
  geom_text(
    data = band_labels_DF10 %>% mutate(y = ifelse(line == 3, 2.1, y)),
    aes(x = Inf, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 1.1,
    fontface = "bold",
    size = 10,
    parse = TRUE
  ) +
  geom_point(
    aes(fill = date, shape = treatment),
    position = pos3,
    size = 6,
    colour = "black",
    stroke = 1.2
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  scale_fill_manual(
    values = date_cols,
    breaks = date_levels,
    drop = FALSE,
    name = NULL
  ) +
  scale_shape_manual(
    values = treatment_shapes,
    breaks = treatment_levels,
    drop = FALSE,
    name = NULL
  ) +
  
  scale_y_log10(
    breaks = c(1, 3, 10, 100),
    limits = c(NA, 320)
  ) +
  scale_x_discrete(
    labels = c(
      "PXR_DF10"   = expression(paste("PXR"[DF10])),
      "ERα_DF10"   = expression(paste("ERα"[DF10])),
      "Nrf2_DF10"   = expression(paste("Nrf2"[DF10])),
      "PFAS_DF10"   = expression(paste("PFAS"[DF10])),
      "PAH_DF10"   = expression(paste("PAH"[DF10]))
    )) +
  
  labs(
    x = NULL,
    y = "EBT-Exceedance"
  ) +
  guides(
    fill = guide_legend(
      title = NULL,
      override.aes = list(
        shape = 21,
        colour = "black",
        size = 5,
        stroke = 1.1
      )
    ),
    shape = guide_legend(
      title = NULL,
      override.aes = list(
        fill = "grey90",
        colour = "black",
        size = 5,
        stroke = 1.1
      )
    )
  ) +
  
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
p1_point_EBT_SW_date

p2_point_EBT_DW_date <- df_EBT_DF_plot %>%
  mutate(calux_assay = factor(calux_assay, levels = assay_levels_DF)) %>% 
  filter(!is.na(values), str_detect(variable, "_dw_")) %>%
  ggplot(aes(x = calux_assay, y = values)
  ) +
  geom_hline(
    yintercept = c(1, 3, 10, 100),
    colour = "grey50",
    linetype = "dashed"
  ) +
  geom_text(
    data = band_labels_DF10 %>% mutate(y = ifelse(line == 1, 0.25, y)),
    aes(x = Inf, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 1.1,
    fontface = "bold",
    size = 10,
    parse = TRUE
  ) +
  geom_point(
    aes(fill = date, shape = treatment),
    position = pos3,
    size = 6,
    colour = "black",
    stroke = 1.2
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  scale_fill_manual(
    values = date_cols,
    breaks = date_levels,
    drop = FALSE,
    name = NULL
  ) +
  scale_shape_manual(
    values = treatment_shapes,
    breaks = treatment_levels,
    drop = FALSE,
    name = NULL
  ) +
  
  scale_y_log10(
    breaks = c(1, 3, 10, 100),
    limits = c(NA, 320)
  ) +
  scale_x_discrete(
    labels = c(
      "PXR_DF10"   = expression(paste("PXR"[DF10])),
      "ERα_DF10"   = expression(paste("ERα"[DF10])),
      "Nrf2_DF10"   = expression(paste("Nrf2"[DF10])),
      "PFAS_DF10"   = expression(paste("PFAS"[DF10])),
      "PAH_DF10"   = expression(paste("PAH"[DF10]))
    )) +
  
  labs(
    x = NULL,
    y = "EBT-Exceedance"
  ) +
  guides(
    fill = guide_legend(
      title = NULL,
      override.aes = list(
        shape = 21,
        colour = "black",
        size = 5,
        stroke = 1.1
      )
    ),
    shape = guide_legend(
      title = NULL,
      override.aes = list(
        fill = "grey90",
        colour = "black",
        size = 5,
        stroke = 1.1
      )
    )
  ) +
  
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
p2_point_EBT_DW_date


## ============================================= #
## DF Sensitivity Analysis                    ####
## ============================================= #

## ============================================= #
## SW                                            #
## ============================================= #
df_range <- seq(1, 150, by = 0.1)
df_sens_SW <- df_EBM_mean_eff_blank %>%
  filter(!is.na(ebt_sw_ngl)) %>%
  select(treatment, date, calux_assay,
         mean_eff, mean_loq_eff, mean_loq2_eff,
         mean_eff_blank_cor,
         mean_loq_eff_blank_cor,
         mean_loq2_eff_blank_cor,
         ebt_sw_ngl)%>%
  pivot_longer(
    cols = c(
      mean_eff,
      mean_loq_eff,
      mean_loq2_eff,
      mean_eff_blank_cor,
      mean_loq_eff_blank_cor,
      mean_loq2_eff_blank_cor
    ),
    names_to = "variable",
    values_to = "value",
    values_drop_na = TRUE
  )

df_sens_SW_range <- df_sens_SW %>%
  crossing(DF = df_range) %>%
  mutate(
    result = value / DF,
    no_exceedance = result <= ebt_sw_ngl
  )

required_DF_SW <- df_sens_SW_range %>%
  group_by(treatment, calux_assay, variable, DF) %>%
  summarise(
    all_results_below_limit = all(no_exceedance),
    maximum_result = max(result),
    .groups = "drop"
  ) %>%
  filter(all_results_below_limit) %>%
  group_by(treatment, calux_assay, variable) %>%
  slice_min(DF, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(
    required_DF = DF,
    result_at_required_DF = maximum_result
  ) %>%
  mutate(matrix = "SW")

required_DF_SW %>% filter(treatment == "GAC") %>%
  slice_max(required_DF, n = 1, with_ties = TRUE)
required_DF_SW %>%
  select(-result_at_required_DF) %>% 
  filter(treatment == "GAC") %>%
  group_by(calux_assay) %>%
  slice_max(required_DF, n = 1, with_ties = TRUE) %>%
  ungroup()

## ============================================= #
## DW                                            #
## ============================================= #
df_sens_DW <- df_EBM_mean_eff_blank %>%
  filter(!is.na(ebt_dw_ngl)) %>%
  select(
    treatment, date, calux_assay,
    mean_eff,
    mean_loq_eff,
    mean_loq2_eff,
    mean_eff_blank_cor,
    mean_loq_eff_blank_cor,
    mean_loq2_eff_blank_cor,
    ebt_dw_ngl
  ) %>%
  pivot_longer(
    cols = c(
      mean_eff,
      mean_loq_eff,
      mean_loq2_eff,
      mean_eff_blank_cor,
      mean_loq_eff_blank_cor,
      mean_loq2_eff_blank_cor
    ),
    names_to = "variable",
    values_to = "value",
    values_drop_na = TRUE
  )

df_sens_DW_range <- df_sens_DW %>%
  crossing(DF = df_range) %>%
  mutate(
    result = value / DF,
    no_exceedance = result <= ebt_dw_ngl
  )

required_DF_DW <- df_sens_DW_range %>%
  group_by(treatment, calux_assay, variable, DF) %>%
  summarise(
    all_results_below_limit = all(no_exceedance),
    maximum_result = max(result),
    .groups = "drop"
  ) %>%
  filter(all_results_below_limit) %>%
  group_by(treatment, calux_assay, variable) %>%
  slice_min(DF, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(
    required_DF = DF,
    result_at_required_DF = maximum_result
  ) %>%
  mutate(matrix = "DW")

required_DF_DW %>%
  select(-result_at_required_DF) %>% 
  filter(treatment == "GAC") %>%
  group_by(calux_assay) %>%
  slice_max(required_DF, n = 1, with_ties = TRUE) %>%
  ungroup()
