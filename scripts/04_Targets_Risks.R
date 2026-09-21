library(tidyverse)
library(janitor)

# ============================================= #
# Read in data                              ####
# ============================================= #

# Read in raw data files
raw_norman <- read_excel("raw_data/risk_analysis/NORMAN_PNEC_20260322.xlsx")
name_lookup <- read_excel("raw_data/risk_analysis/replace_names.xlsx")
raw_data_pgv <- read_excel("raw_data/risk_analysis/pgv_targets_2025-10-18.xlsx")

# =============================================== #
# Data Wrangling                               ####
# =============================================== #

## ============================================== #
## NORMAN data                                 ####
## ============================================== #

# Clean names, convert pnec to ng/l, rename columns, rename analytes, transform µg/L to ng/L
df_norman <- raw_norman %>%
  rename(compound = Compound) %>%
  mutate(
    compound = case_when(
      compound %in% name_lookup$old_name ~ name_lookup$new_name[match(compound, name_lookup$old_name)],
      TRUE ~ compound
    ),
    pnec_ngL = as.numeric(pnec_ugL) * 1000
  ) %>%
  select(-pnec_ugL)

## ============================================== #
## pGv data                                    ####
## ============================================== #
# Note: pGV = (preliminary) guideline values for drinking water


# Select columns, rename value column, clean names, rename analytes
df_pgv <- raw_data_pgv %>%
  select(c(Compound, CID, CanonicalSMILES, value, unit)) %>%
  rename(pgv_ngL = value) %>%
  select(-c("CID", "CanonicalSMILES", "unit")) %>%
  clean_names()
df_pgv <- df_pgv %>% mutate(compound = case_when(
  compound %in% name_lookup$old_name ~ name_lookup$new_name[match(compound, name_lookup$old_name)],
  TRUE ~ compound
))

## ============================================== #
## OMP target data                             ####
## ============================================== #

# Assign to new df, convert class to title case
df_OMP_risks <- df_OMP %>% # df_OMP was generated in the script 01_Target_Analysis
  select(c(analyte_name, sample_date, treatment, replicate, quantity_EF, class)) %>% mutate(class = str_to_title(class)) %>% 
  mutate(analyte_name = ifelse(analyte_name == "1-Phenylurea", "N-Phenylurea", analyte_name))

# ============================================= #
# Data Analysis                              ####
# ============================================= #

## ============================================= #
## Env. Risks, Freshwater                     ####
## ============================================= #


# Join OMP risks data with NORMAN PNEC data
df_risks_env <- full_join(df_OMP_risks, df_norman, by = c("analyte_name" = "compound"))

# Calculate mean RQs per day and treatment
df_env_risks_mean <- df_risks_env %>%
  group_by(analyte_name, sample_date, treatment, replicate) %>%
  mutate(
    RQ = quantity_EF / pnec_ngL,
    RQ_DF10 = RQ / 10
  ) %>%
  ungroup() %>%
  group_by(analyte_name, sample_date, treatment) %>%
  summarise(
    mean_RQ = mean(RQ, na.rm = TRUE),
    sd_RQ = sd(RQ, na.rm = TRUE),
    rsd_RQ = sd_RQ / mean_RQ * 100,
    mean_RQ_DF10 = mean(RQ_DF10, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(desc(mean_RQ))

# Calculate total RQs per treatment and day
df_risk_removal <- df_env_risks_mean %>% 
  filter(!is.na(sample_date)) %>% 
  group_by(sample_date, treatment) %>% 
  summarise(sum_RQ = sum(mean_RQ, na.rm = TRUE)) %>% 
  group_by(sample_date) %>% 
  mutate(baseline = sum_RQ[treatment=="WWTP-E"],
         risk_rem = (sum_RQ - baseline)/baseline * -100) %>% 
  select(-baseline)
df_risk_removal %>% filter(treatment == "GAC" | sample_date == "19.07.2022" & treatment == "CMF")

# Plot total RQs per treatment and day
total_RQ <- df_risk_removal %>% 
  ggplot(aes(x = treatment, y = sum_RQ)) +
  geom_point() +
  scale_y_log10() +
  labs(x = NULL,
       y = "Total RQ") +
  facet_wrap(. ~ sample_date) +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold"),
    axis.text.x = element_text(angle = 40, hjust = 1)
  )
total_RQ


# Which compounds are above RQ 0.1 after GAC?
df_env_risks_mean %>%
  filter(mean_RQ >= 0.1) %>%
  filter(treatment == "GAC") %>%
  distinct(analyte_name) %>%
  arrange(desc(.))

# On how many days are those compounds above RQ 0.1 after GAC?
df_env_risks_mean %>%
  filter(mean_RQ >= 0.1) %>%
  filter(treatment == "GAC") %>%
  group_by(analyte_name, treatment) %>%
  summarise(n = length(sample_date)) %>%
  arrange(desc(n))

# Filter for OMPs above RQ >= 1 & count number of exceeded days
df_env_risks_mean %>%
  filter(mean_RQ >= 1) %>%
  group_by(analyte_name, treatment) %>%
  summarise(n = length(sample_date)) %>%
  arrange(desc(n))

# Prepare data for supplementary information (SW)
tabSI_SW <- left_join(df_env_risks_mean, select(df_norman, -2, -3), by = c("analyte_name" = "compound"), relationship = "many-to-many") %>%
  distinct(analyte_name, sample_date, treatment, .keep_all = TRUE) %>%
  filter(mean_RQ != 0) %>% 
  relocate(pnec_ngL, .after = mean_RQ_DF10)

# Calculate mean RQ for compound x per treatment
comp_x <- "Venlafaxine"
df_env_risks_mean %>%
  filter(analyte_name == comp_x) %>%
  group_by(treatment) %>%
  summarise(
    mean_RQ = mean(mean_RQ, na.rm = TRUE),
    mean_RQ_DF10 = mean(mean_RQ_DF10, na.rm = TRUE)
  )


## ============================================= #
## Drinking Water Risks                       ####
## ============================================= #

# Join OMP risks data with pGV data
df_risks_dw <- left_join(df_OMP_risks, df_pgv, by = c("analyte_name" = "compound"))

# Number of compounds with pgv
df_risks_dw %>%
  distinct(analyte_name, .keep_all = TRUE) %>%
  reframe(
    n_analytes = length(analyte_name),
    n_pgv = sum(!is.na(pgv_ng_l)),
    n_NA = sum(is.na(pgv_ng_l)),
    frac_pgv = round(n_pgv / n_analytes * 100, digits = 0)
  )

# Calculate mean RQs
df_dw_RQ_mean <- df_risks_dw %>%
  group_by(analyte_name, sample_date, treatment) %>%
  mutate(
    RQ = quantity_EF / pgv_ng_l,
    RQ_DF10 = RQ / 10
  ) %>%
  ungroup() %>%
  group_by(analyte_name, sample_date, treatment) %>%
  summarise(
    mean_RQ = round(mean(RQ, na.rm = TRUE), digits = 4),
    sd_RQ = round(sd(RQ, na.rm = TRUE), digits = 4),
    rsd_RQ = round(sd_RQ / mean_RQ * 100, digits = 4),
    mean_RQ_DF10 = round(mean(RQ_DF10, na.rm = TRUE), digits = 4)
  ) %>%
  ungroup() %>%
  arrange(desc(mean_RQ))

# Filter for OMPs above RQ >= 1 & count number of exceeded days
df_dw_RQ_mean %>%
  filter(mean_RQ >= 1) %>%
  group_by(analyte_name, treatment) %>%
  summarise(n = length(sample_date)) %>%
  arrange(desc(n))


# Prepare data for supplementary information (DW)
tabSI_DW <- left_join(df_dw_RQ_mean, select(df_risks_dw, 1, 7), by = "analyte_name", relationship = "many-to-many") %>%
  distinct(analyte_name, sample_date, treatment, .keep_all = TRUE) %>%
  filter(mean_RQ != 0)


## ============================================= #
## Source water for drinking water            ####
## ============================================= #

df_source_w <- df_sum_analyte %>% # df_sum_analyte was generated in the script 01_Target_Analysis
  mutate(
    SW_limit_ngl = 1000,
    mean_SW_exceedance = mean / SW_limit_ngl
  ) %>%
  arrange(desc(mean_SW_exceedance)) %>%
  select(c(sample_date, treatment, analyte_name, mean_SW_exceedance))

# How many compounds per day exceed the threshold
df_source_w %>%
  group_by(sample_date, treatment) %>%
  summarise(n_above1 = sum(mean_SW_exceedance > 1)) %>%
  pivot_wider(names_from = treatment, values_from = n_above1)

# Which compound exceeds the threshold on day_x
day_x <- "19.07.2022"
df_source_w %>% filter(mean_SW_exceedance > 1, sample_date == day_x)

# Which compounds exceed the threshold in the wastewater treatment plant effluent
df_source_w %>%
  filter(mean_SW_exceedance > 1, treatment == "WWTP-E") %>%
  distinct(analyte_name)