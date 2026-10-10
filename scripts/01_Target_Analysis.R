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
                        "26.08.2022", "15.09.2022", "19.10.2022")
treatment_order <- c("WWTP-E", "O3", "AO", "CMF", "GAC")
df_OMP_raw <- df_OMP_raw %>% mutate(
  sample_date = factor(sample_date, levels = sample_dates_order),
  treatment = factor(treatment, levels = treatment_order)
)

# Filter for detected analytes and samples
df_OMP <- df_OMP_raw %>%
  filter(TASQ_sample_type == "Sample") %>%
  filter(!is.na(sample_date)) 

# Which OMP classes were detected?
df_OMP %>% 
  filter(quantity_raw != 0, TASQ_sample_type == "Sample") %>% 
  group_by(class) %>%
  distinct(analyte_name, .keep_all = TRUE) %>%
  summarise(n = n()) %>%
  mutate(total = sum(n),
         contri = n / total * 100)

# ============================================= #
# Summary statistics                         #### 
# ============================================= #

# Summary per treatment, date, and analyte
df_sum_analyte <- df_OMP %>%
  mutate(quant.for.summary = ifelse(limit != "quantity", 0, quantity_EF)) %>% 
  filter(limit != ">C1") %>% 
  group_by(sample_date, treatment, analyte_name) %>%
  summarise(
    n = n(),
    n_above0 = sum(quant.for.summary > 0),
    min = min(quant.for.summary, na.rm = TRUE),
    max = max(quant.for.summary, na.rm = TRUE),
    median = median(quant.for.summary, na.rm = TRUE),
    iqr = IQR(quant.for.summary, na.rm = TRUE),
    mean = mean(quant.for.summary, na.rm = TRUE),
    sd = sd(quant.for.summary, na.rm = TRUE),
    rsd = sd / mean * 100,
    .groups = "drop"
  )

# Number of compounds per day and treatment that could be quantified
df_n_comp_conc <- df_sum_analyte %>% group_by(treatment, sample_date) %>% 
  summarise(n_quant = sum(mean > 0, na.rm = TRUE)) %>% 
  ungroup()

# ============================================= #
# Removal Calculation                        ####
# ============================================= #

# Total mean concentration of all analytes
df_rem_total <- df_OMP %>%
  mutate(quant.for.summary = ifelse(limit != "quantity", 0, quantity_EF)) %>% 
  filter(limit != ">C1") %>% 
  group_by(sample_date, treatment, analyte_name) %>%
  summarise(
    mean = mean(quant.for.summary, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
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

# Removal per treatment, date, and analyte
df_rem_analyte_day <- df_OMP %>%
  distinct(treatment, sample_date, analyte_name) %>%
  left_join(
    df_sum_analyte %>% select(c(1:3, mean)),
    by = c("treatment", "sample_date", "analyte_name")
  ) %>%
  arrange(sample_date, analyte_name, treatment) %>%
  group_by(sample_date, analyte_name) %>%
  mutate(
    baseline = mean[treatment == "WWTP-E"][1],
    rel_removal_overall = (1 - mean / baseline) * 100,
    rel_removal_overall = ifelse(rel_removal_overall == 0, NA, rel_removal_overall),
    rel_removal_overall = round(rel_removal_overall),
    rel_removal_overall = ifelse(is.infinite(rel_removal_overall), -100, rel_removal_overall)
  ) %>%
  arrange(treatment, .by_group = TRUE) %>%
  mutate(
    prev_conc = lag(mean),
    rel_removal_treatment = (1 - mean / prev_conc) * 100
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


# ============================================= #
# Trend Classification                       ####
# ============================================= #

# Assign treatment trends
df_rem_trend <- df_rem_analyte_day %>%
  group_by(sample_date, analyte_name) %>%
  arrange(treatment, .by_group = TRUE) %>%
  mutate(
    trend = case_when(
      treatment == "WWTP-E" ~ "WWTP-E",
      
      is.na(mean) & is.na(prev_conc) ~ "Complete data missing (>C1)",
      is.na(mean) ~ "Complete data missing (>C1)",
      is.na(prev_conc) ~ "Complete data missing (>C1)",
      
      mean == 0 & prev_conc == 0 & baseline > 0 ~ "Removed Previously",
      mean == 0 & prev_conc == 0 & baseline == 0 ~ "Not detected",
      mean == 0 & prev_conc == 0 ~ "Complete data missing (>C1)",
      
      rel_removal_treatment > 0 &
        rel_removal_treatment < 100 ~ "Reduced",
      
      rel_removal_treatment < 0 &
        is.finite(rel_removal_treatment) ~ "Increased",
      
      rel_removal_treatment == -Inf ~ "Newly detected",
      rel_removal_treatment == 100 ~ "Removed",
      rel_removal_treatment == 0 ~ "Unchanged",
      
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()


# Count removal trends per sampling date and treatment
df_n_comp_rem <- df_rem_trend %>%
  group_by(sample_date, treatment) %>%
  summarise(
    n_reduced = sum(trend == "Reduced", na.rm = TRUE),
    n_increased = sum(trend == "Increased", na.rm = TRUE),
    n_new_detected = sum(trend == "Newly detected", na.rm = TRUE),
    n_removed = sum(trend == "Removed", na.rm = TRUE),
    n_unchanged = sum(trend == "Unchanged", na.rm = TRUE),
    n_still_removed = sum(trend == "Removed Previously", na.rm = TRUE),
    n_not_detected = sum(trend == "Not detected", na.rm = TRUE),
    n.current.unavailable = sum(trend == "Complete data missing (>C1)", na.rm = TRUE),
    n.previous.unavailable = sum(trend == "Complete data missing (>C1)", na.rm = TRUE),
    n.both.unavailable = sum(trend == "Complete data missing (>C1)", na.rm = TRUE),
    n.zero.baseline.unavailable = sum(trend == "Complete data missing (>C1)", na.rm = TRUE),

    # >C1 compounds. Do not add to total count, only diagnostic
    n.C1 = sum(is.na(rel_removal_treatment) & is.na(baseline) & !is.na(mean), na.rm = TRUE),

    # Positive concentrations with an unavailable comparison step
    n_cur_pos_without_prev_C1 = sum(trend == "Complete data missing (>C1)" & mean > 0, na.rm = TRUE),
    n_prev_pos_without_cur_C1 = sum(trend == "Complete data missing (>C1)" & prev_conc > 0, na.rm = TRUE),

    # Unique compounds with at least one available removal metric.
    n_overall_removal_val = n_distinct(
      analyte_name[!is.na(rel_removal_overall) | !is.na(rel_removal_treatment)],
      na.rm = TRUE
    ),

    # Analytes accounted for by trends and availability flags
    n.current.accounted =
      n_reduced + n_increased + n_new_detected + n_unchanged + n_cur_pos_without_prev_C1,
    n.previous.accounted = n_reduced + n_increased + n_removed + n_unchanged + n_prev_pos_without_cur_C1,
    .groups = "drop"
  )

# Combine concentration/removal counts and compare
df_n_comp <- df_n_comp_conc %>%
  left_join(df_n_comp_rem, by = c("treatment", "sample_date")) %>%
  group_by(sample_date) %>%
  arrange(treatment, .by_group = TRUE) %>%
  mutate(
    n.previous.quant = lag(n_quant),
    check_current_count = if_else(
      treatment == "WWTP-E",
      NA,
      as.integer(n.current.accounted == n_quant)
    ),
    check_previous_count = if_else(
      treatment == "WWTP-E",
      NA,
      as.integer(n.previous.accounted == n.previous.quant)
    ),
    n_overall_removal_val = if_else(
      treatment == "WWTP-E",
      NA,
      n_overall_removal_val
    )
  ) %>%
  select(-contains(".")) %>%
  ungroup() %>%
  mutate(
    across(
      -all_of(c("treatment", "sample_date", "n_quant")),
      ~ replace(.x, treatment == "WWTP-E", NA)
    )
  )

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
  filter(fod_detection >= 1) %>%
  print(n = nrow(.))

# Which compounds were found with at least X % detection frequency for individual treatments?
detection_treatment %>%
  filter(fod_detection >= 1) %>%
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

df_n_comp2 <-  df_sum_analyte %>%
  group_by(treatment, sample_date) %>%
  summarise(max_val = max(mean, na.rm = TRUE),
            .groups = "drop"
  ) %>%
  mutate(max_val = if_else(treatment == "O3" & sample_date == "26.07.2022", max_val + 350, max_val)) %>% 
  right_join(., df_n_comp, by = c("sample_date", "treatment"))

pos <- position_dodge(width = 0.8)

p1_SI <- ggplot(
  data = df_sum_analyte,
  aes(x = treatment, y = mean, fill = sample_date)) +
  geom_boxplot(position = pos) +
  geom_text(data = df_n_comp2,
            aes(x = treatment, y = max_val, label = n_quant),
            vjust = -0.5,
            position = pos,
            size = 4
  ) +
  labs(
    x = NULL,
    y = expression(paste("Concentration (ngL"^-1, ")")),
    fill = NULL
  ) +
  scale_y_log10(labels = label_log(digits = 2)) +
  coord_cartesian(ylim = c(-Inf, max(df_n_comp2$max_val) + 1000)) +
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
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.key.size = unit(1, "cm"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey98")
  )
p1_SI

# =========================================================== #
# Boxplot of removal efficiency per treatment and date

p2_SI <- df_rem_analyte_day %>%
  filter(treatment != "WWTP-E") %>%
  ggplot(aes(x = treatment, y = rel_removal_overall, fill = sample_date)) +
  geom_boxplot(position = pos) +
  geom_text(data = filter(df_n_comp, treatment != "WWTP-E"),
            aes(x = treatment, y = 107, label = n_overall_removal_val),
            position = pos,
            size = 4
  ) +
  scale_fill_brewer(palette = "RdBu") +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  coord_cartesian(ylim = c(-100, 110)) +
  labs(
    x = NULL,
    y = "Removal efficiency (%)",
    fill = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.key.size = unit(1, "cm"),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "pt"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey98")
  )
p2_SI

# =========================================================== #
# Scatterplot of relative removal per analyte and treatment

# Order based on GAC total removal
order_analytes_pre <- df_rem_analyte_total %>%
  filter(treatment == "GAC") %>%
  mutate(rel_removal_overall = ifelse(is.na(rel_removal_overall), 0, rel_removal_overall)) %>% 
  arrange(rel_removal_overall) %>%
  pull(analyte_name)

# Compounds without GAC removal data will end up as NAs in the final plot. Combine with all names again
all_analytes <- df_rem_analyte_total %>%
  distinct(analyte_name) %>%
  pull(analyte_name)
order_analytes <- c(setdiff(order_analytes_pre, all_analytes), order_analytes_pre)

figure_S3 <- df_rem_analyte_total %>%
  filter(treatment != "WWTP-E", !is.na(rel_removal_overall)) %>%
  mutate(analyte_name = factor(analyte_name, levels = order_analytes)) %>%
  ggplot(aes(y = analyte_name, x = rel_removal_overall, col = treatment)) +
  geom_vline(aes(xintercept = 80), linetype = "dashed") +
  geom_point(alpha = 0.8) +
  coord_cartesian(xlim = c(-50, 110)) +
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

# =========================================================== #
# Cumulative concentrations per class - stacked barplot 
figure_S19 <- df_OMP %>%
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
figure_S19