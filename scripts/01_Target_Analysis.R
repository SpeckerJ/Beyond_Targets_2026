library(readxl)
library(tidyverse)
library(scales)

# ============================================= #
# Read in data                               ####
# ============================================= #

df_OMP_raw <- read_excel("raw_data/OMP/raw_OMP.xlsx")



# ============================================= #
# Data wrangling                             ####
# ============================================= #

# Change to quantity_raw 0 for ND/LOQ/LOD
df_OMP_raw <- df_OMP_raw %>% mutate(
  quantity_raw = case_when(
    limit == "ND" ~ 0,
    limit == "<LOD" ~ 0,
    limit == "<LOQ" ~ 0,
    limit == "quantity" ~ quantity_raw,
    limit == ">C1" ~ quantity_raw,
    TRUE ~ quantity_raw
  )
)

# Blank subtraction
df_OMP_raw <- df_OMP_raw %>%
  group_by(analyte_name) %>%
  mutate(
    mean_conc_blank = mean(quantity_raw[TASQ_sample_type == "Blank"], na.rm = TRUE),
    quantity_blank_corr = ifelse(quantity_raw > 0 & TASQ_sample_type == "Sample",
      quantity_raw - mean_conc_blank, quantity_raw
    )
  ) %>%
  ungroup()

# Enrichment Factor
df_OMP_raw <- df_OMP_raw %>%
  group_by(analyte_name) %>%
  mutate(quantity_EF = quantity_blank_corr / EF) %>%
  ungroup()

# Order sample dates and treatments for consistent plotting
sample_dates_order <- c("19.07.2022", "26.07.2022", "16.08.2022",
                        "26.08.2022", "15.09.2022", "15.10.2022")
treatment_order <- c("WWTP-E", "O3", "AO", "CMF", "GAC")
df_OMP_raw <- df_OMP_raw %>% mutate(
  sample_date = factor(sample_date, levels = sample_dates_order),
  treatment = factor(treatment, levels = treatment_order)
)

# Filter for detected analytes and samples
df_OMP <- df_OMP_raw %>%
  filter(TASQ_sample_type == "Sample") %>%
  filter(!is.na(sample_date))

# ============================================= #
# Summary statistics                         #### 
# ============================================= #

# Summary per treatment and date
df_sum_overall <- df_OMP %>%
  filter(limit == "quantity") %>% # Filter for samples within the calibration line; exclude ND, <LOD, <LOQ
  group_by(sample_date, treatment) %>%
  summarise(
    n = sum(quantity_EF > 0),
    min_conc = min(quantity_EF, na.rm = TRUE),
    max_conc = max(quantity_EF, na.rm = TRUE),
    median_conc = median(quantity_EF, na.rm = TRUE),
    iqr_conc = IQR(quantity_EF, na.rm = TRUE),
    mean_conc = mean(quantity_EF, na.rm = TRUE),
    sd_conc = sd(quantity_EF, na.rm = TRUE),
    rsd_conc = sd_conc / mean_conc * 100,
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(., digits = 2)))


# Summary per treatment, date, and analyte
df_sum_analyte <- df_OMP %>%
  filter(limit == "quantity") %>%
  group_by(sample_date, treatment, analyte_name) %>%
  summarise(
    n = sum(quantity_EF > 0),
    min = min(quantity_EF, na.rm = TRUE),
    max = max(quantity_EF, na.rm = TRUE),
    median = median(quantity_EF, na.rm = TRUE),
    iqr = IQR(quantity_EF, na.rm = TRUE),
    mean = mean(quantity_EF, na.rm = TRUE),
    sd = sd(quantity_EF, na.rm = TRUE),
    rsd = sd / mean * 100,
    .groups = "drop"
  )

# ============================================= #
# Removal Calculation                        ####
# ============================================= #

# Total mean concentration of all analytes
df_rem_total <- df_sum_analyte %>%
  group_by(sample_date, treatment) %>%
  summarise(sum_conc = round(sum(mean))) %>%
  group_by(sample_date) %>%
  mutate(
    baseline = sum_conc[treatment == "WWTP-E"],
    rel_removal_total = round((1 - sum_conc / baseline) * 100),
    rel_removal_total = ifelse(rel_removal_total == 0, NA, rel_removal_total)
  ) %>%
  ungroup() %>%
  select(-baseline)


# Removal per treatment and date
df_rem <- df_OMP %>%
  filter(limit == "quantity") %>%
  group_by(treatment, sample_date) %>%
  summarise(
    mean_conc = mean(quantity_EF, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(sample_date) %>%
  mutate(
    baseline = mean_conc[treatment == "WWTP-E"],
    rel_removal_overall = (1 - mean_conc / baseline) * 100,
    rel_removal_overall = ifelse(rel_removal_overall == 0, NA, rel_removal_overall),
    rel_removal_overall = round(rel_removal_overall)
  ) %>%
  ungroup()

# Removal per treatment, date, and analyte
df_rem_analyte_day <- df_OMP %>%
  group_by(treatment, sample_date, analyte_name) %>%
  summarise(
    mean_conc = round(mean(quantity_EF, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  group_by(sample_date, analyte_name) %>%
  mutate(
    baseline = mean_conc[treatment == "WWTP-E"][1],
    rel_removal_overall = (1 - mean_conc / baseline) * 100,
    rel_removal_overall = ifelse(rel_removal_overall == 0, NA, rel_removal_overall),
    rel_removal_overall = round(rel_removal_overall)
  ) %>%
  arrange(treatment, .by_group = TRUE) %>%
  mutate(
    prev_conc = lag(mean_conc),
    rel_removal_treatment = (1 - mean_conc / prev_conc) * 100
  ) %>%
  ungroup()

# Total removal per analyte across all dates
df_rem_analyte_total <- df_OMP %>%
  group_by(treatment, analyte_name) %>%
  summarise(
    mean_conc = round(mean(quantity_EF, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  group_by(analyte_name) %>%
  mutate(
    baseline = mean_conc[treatment == "WWTP-E"][1],
    rel_removal_overall = (1 - mean_conc / baseline) * 100,
    rel_removal_overall = ifelse(rel_removal_overall == 0, NA, rel_removal_overall),
    rel_removal_overall = round(rel_removal_overall),
    rel_removal_overall = ifelse(is.infinite(rel_removal_overall) & rel_removal_overall < 0, -100, rel_removal_overall)
  ) %>%
  arrange(treatment, .by_group = TRUE) %>%
  mutate(
    prev_conc = lag(mean_conc),
    rel_removal_treatment = (1 - mean_conc / prev_conc) * 100
  ) %>%
  ungroup()


# Analyte removal rates #
df_rem_80 <- df_rem_analyte_total %>% filter(rel_removal_overall >= 80)
df_rem_20 <- df_rem_analyte_total %>% filter(rel_removal_overall <= 20)

# ============================================ #
# Frequency of Detection                    ####
# ============================================ #

# Define detection/non-detection
analyte_detected <- c("quantity", "<LOQ", ">C1")
analyte_nondetect <- c("<LOD", "ND")

# Calculate the detection/non-detection for each analyte per date and treatment
limit_abs_det1 <- df_OMP %>%
  group_by(analyte_name, sample_date, treatment, limit) %>%
  summarise(n_limit = n(), .groups = "drop")

# Fill missing combinations and add 'detected' flag
limit_abs_det2 <- limit_abs_det1 %>%
  complete(analyte_name, nesting(sample_date, treatment, limit), fill = list(n_limit = 0)) %>%
  mutate(
    detected = case_when(
      limit %in% analyte_detected ~ "1",
      limit %in% analyte_nondetect ~ "0",
      TRUE ~ NA_character_
    ),
    detected = as.integer(detected)
  )

# Calculate the overall detection frequency over the entire sampling campaign
detection_overall <- limit_abs_det2 %>%
  group_by(analyte_name, sample_date) %>%
  filter(detected == 1, n_limit >= 1) %>%
  summarise(
    detected_day = as.integer(any(detected == 1)), .groups = "drop"
  ) %>%
  complete(analyte_name, nesting(sample_date), fill = list(detected_day = 0)) %>%
  group_by(analyte_name) %>%
  summarise(
    fod_detection = sum(detected_day) / length(sample_date),
    fod_nondect = 1 - fod_detection,
    .groups = "drop"
  )
detection_overall

# Calculate the detection frequency over the entire sampling campaign per treatment
detection_treatment <- limit_abs_det2 %>%
  group_by(analyte_name, sample_date, treatment) %>%
  filter(detected == 1, n_limit >= 1) %>%
  summarise(
    detected_day = as.integer(any(detected == 1)), .groups = "drop"
  ) %>%
  complete(analyte_name, nesting(sample_date, treatment), fill = list(detected_day = 0)) %>%
  group_by(analyte_name, treatment) %>%
  summarise(
    fod_detection = sum(detected_day) / length(sample_date),
    fod_nondect = 1 - fod_detection,
    .groups = "drop"
  )
detection_treatment

# Which compounds were found with at least X % detection frequency over the entire sampling campaign/on every sampling day?
detection_overall %>%
  filter(fod_detection == 1) %>%
  print(n = nrow(.))

# Which compounds were found with at least X % detection frequency for individual treatments?
detection_treatment %>%
  filter(fod_detection == 1) %>%
  print(n = nrow(.))

# Which compounds were always found over the entire treatment train/in every sample?
detection_treatment %>%
  filter(fod_detection == 1) %>%
  group_by(analyte_name) %>%
  filter(length(treatment) == 5) %>%
  pull(analyte_name) %>%
  unique()

# How many compounds were always found over the entire treatment train/in every sample?
detection_treatment %>%
  filter(fod_detection == 1) %>%
  group_by(analyte_name) %>%
  filter(length(treatment) == 5) %>%
  pull(analyte_name) %>%
  unique() %>%
  length()


# ============================================= #
# Plots                                      ####
# ============================================= #

# Set theme for all (!) following plots
theme_set(
  theme_bw(base_size = 26) +
    theme(
      # Text color
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      
      # Panel border
      panel.border = element_rect(color = "black", fill = NA),
      strip.background = element_rect(fill = "white"),
      strip.text.x = element_text(face = "bold")
    )
)

# =========================================================== #
# Boxplot mean concentration per treatment and date

pos <- position_dodge(width = 0.8)

p1_SI <- ggplot(
  data = df_sum_analyte,
  aes(x = treatment, y = mean, fill = sample_date)
) +

  # boxplots with dodging
  geom_boxplot(position = pos) +
  geom_point(
    data = df_rem_total,
    aes(x = treatment, y = sum_conc, fill = sample_date, group = sample_date),
    position = pos,
    shape = 22,
    size = 5,
    show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = expression(paste("Concentration (ngL"^-1, ")")),
    fill = NULL
  ) +
  scale_y_log10(labels = label_log(digits = 2)) +
  scale_fill_brewer(palette = "RdBu") +
  scale_x_discrete(
    drop = FALSE,
    labels = c(
      "O3" = expression(paste("O"[3]))
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    legend.position = "bottom"
  )
p1_SI

# =========================================================== #
# Boxplot of removal efficiency per treatment and date

p2_SI <- df_rem_analyte_day %>%
  filter(treatment != "WWTP-E") %>%
  ggplot(aes(x = treatment, y = rel_removal_overall, fill = sample_date)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "RdBu") +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  coord_cartesian(ylim = c(-100, 100)) +
  labs(
    x = NULL,
    y = "Removal efficiency (%)",
    fill = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    legend.position = "bottom"
  )
p2_SI

# =========================================================== #
# Cumulative concentrations per class - stacked barplot 

OMP_pesticides <- df_OMP %>%
  filter(class %in% c("pesticide")) %>%
  mutate(fill_arg = case_when(
    analyte_name %in% c("DEET", "MCPA") ~ analyte_name,
    TRUE ~ "Other pesticides"
  )) %>%
  ggplot(aes(x = treatment, y = quantity_EF, fill = fill_arg)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_hline(yintercept = 500, col = "black", linetype = "dashed") +
  labs(
    x = NULL,
    y = expression(paste("Concentration (ngL"^-1, ")")),
    fill = NULL,
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  facet_wrap(. ~ sample_date) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(1, "cm"),
    legend.position = "bottom"
  )
OMP_pesticides

# =========================================================== #
# Scatterplot of relative removal per analyte and treatment

# Order based on GAC total removal
order_analytes_pre <- df_rem_analyte_total %>%
  filter(treatment == "GAC") %>%
  arrange(rel_removal_overall) %>%
  pull(analyte_name)

# Compounds without GAC removal data will end up as NAs in the final plot. Combine with all names again
all_analytes <- df_rem_analyte_total %>%
  distinct(analyte_name) %>%
  pull(analyte_name)
order_analytes <- c(setdiff(all_analytes, order_analytes_pre), order_analytes_pre)

figure_S3 <- df_rem_analyte_total %>%
  filter(treatment != "WWTP-E", !is.na(rel_removal_overall)) %>%
  filter(rel_removal_overall > -50) %>%
  mutate(analyte_name = factor(analyte_name, levels = order_analytes)) %>%
  ggplot(aes(y = analyte_name, x = rel_removal_overall, col = treatment)) +
  geom_vline(aes(xintercept = 80), linetype = "dashed") +
  geom_point() +
  scale_color_brewer(
    palette = "Set2",
    labels = c(
      "O3" = expression(paste("O"[3]))
    )
  ) +
  labs(
    x = "Relative Removal (%)",
    y = NULL,
    col = NULL
  )
figure_S3
