library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(forcats)


## ============================================= #
# Main Manuscript                             ####
## ============================================= #

## ============================================== #
## Read in data                                ####
## ============================================== #

# List all files in directory, identify desired file, read in
raw_data_folder <- list.files("raw_data/NTS", full.names = TRUE)
raw_data_file_ID4 <- raw_data_folder[str_detect(raw_data_folder, "ID1234")]
raw_data_file_NTS <- raw_data_folder[str_detect(raw_data_folder, "total")]
raw_data_ID4 <- read.csv(raw_data_file_ID4)
raw_data_NTS <- read.csv(raw_data_file_NTS)

# assign to working df
df_raw_ID4 <- raw_data_ID4
df_raw_NTS <- raw_data_NTS

# Lookup tables; date and treatment order
# Note: The WWTP effluent is named here WWTP_E and not WWTP-E to avoid issues with column names
lookup_dates <- data.frame(
  day = c(paste0("D", 1:6)),
  date = paste0(c("19.07.", "26.07.", "16.08.", "26.08.", "15.09.", "15.10."), "2022")
)
date_order <- lookup_dates$date
lookup_treatmeans <- data.frame(
  s = paste0("S", rep(1:12)),
  treatment = paste0(rep(c(
    "WWTP_E",
    "O3",
    "CMF",
    "GAC"
  ), each = 3))
)
treat_order <- c("WWTP_E", "O3", "AO", "CMF", "GAC")

## ============================================== #
## Data Wrangling                              ####
## ============================================== #

# Assign features measured in both polarities to one polarity
# to avoid unexpected errors during data wrangling; pivot longer
df_nts <- df_raw_NTS %>%
  mutate(
    # mean intensities
    mean_int_pos = rowMeans(across(matches("_Pos_"), ~ replace_na(.x, 0)), na.rm = TRUE),
    mean_int_neg = rowMeans(across(matches("_Neg_"), ~ replace_na(.x, 0)), na.rm = TRUE),

    # assign preferred polarity
    assigned_pol_mean = if_else(mean_int_pos >= mean_int_neg, "Pos", "Neg")
  ) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "sample",
    values_to = "intensity"
  )

# Extract day, treatment and polarity from raw sample name
df_nts <- df_nts %>%
  mutate(
    day = str_extract(sample, "D[1-6]"),
    s = str_extract(sample, "S(?:1[0-2]|[1-9])"),
    polarity = str_extract(sample, "(?<=_)(Pos|Neg)(?=_)")
  ) %>%
  left_join(lookup_dates, by = "day") %>%
  left_join(lookup_treatmeans, by = "s") %>%
  mutate(treatment = case_when(
    day %in% c("D5", "D6") & s %in% c("S4", "S5", "S6") ~ "AO",
    treatment == "Eff" ~ "WWTP_E",
    TRUE ~ treatment
  )) %>%
  select(-day, -s)

# Assign which polarity to keep for features measured in both polarities
# and general data wrangling
df_nts <- df_nts %>%
  left_join(lookup_dates, by = "date") %>%
  filter(polarity == assigned_pol_mean) %>%
  mutate(
    treatment = factor(treatment, levels = treat_order),
    date = factor(date, levels = date_order)
  ) %>%
  select(-neutralMass, -matches("pol_count|pol_sum|sum_|mean_")) %>%
  relocate(c(sample, day, treatment, polarity), .after = group)

# Extract suspect and TP features
suspect_features <- df_raw_ID4 %>%
  select(group) %>%
  distinct() %>%
  mutate(set = "suspect")
TP_features <- df_nts %>%
  filter(str_detect(susp_name, "-TP")) %>%
  select(group) %>%
  distinct() %>%
  mutate(set = "TP")

# Filter for transformation products (TPs)
df_TPs <- df_nts %>%
  group_by(date, group, assigned_pol_mean, treatment) %>%
  mutate(intensity = mean(intensity, na.rm = TRUE)) %>%
  distinct(date, group, assigned_pol_mean, treatment, .keep_all = TRUE) %>%
  select(group, date, treatment, assigned_pol_mean, intensity, susp_name) %>%
  filter(str_detect(susp_name, regex("-[A-Za-z0-9]*TP\\d+", ignore_case = TRUE)))


## ============================================== #
## Trends Wide - Preparation                   ####
## ============================================== #

# Prepare data for wide-format trend analysis
df_prepared <- df_nts %>%
  select(group, ret, mz, date, treatment, assigned_pol_mean, intensity)

# Calculate mean intensity per group, date, treatment, and polarity
df_mean_intensity <- df_prepared %>%
  group_by(date, group, assigned_pol_mean, treatment) %>%
  mutate(intensity = mean(intensity, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(date, group, assigned_pol_mean, treatment, .keep_all = TRUE)

# Identify first and last observation per group and date
df_obs <- df_mean_intensity %>%
  arrange(date, group, assigned_pol_mean, treatment) %>%
  group_by(date, group, assigned_pol_mean) %>%
  mutate(
    present = ifelse(intensity > 0, TRUE, FALSE),
    first_obs = if (any(intensity > 0, na.rm = TRUE)) {
      treatment[intensity > 0][which.min(as.numeric(treatment[intensity > 0]))]
    } else {
      NA
    },
    last_obs = if (any(intensity > 0, na.rm = TRUE)) {
      treatment[intensity > 0][which.max(as.numeric(treatment[intensity > 0]))]
    } else {
      NA
    }
  ) %>%
  ungroup()

# Check last treatment step and presence in WWTP_E or last treatment
df_last_treat <- df_obs %>%
  group_by(date, group, assigned_pol_mean) %>%
  mutate(
    # Check last treatment step. Depends on the day
    last.treat = last(treatment),

    # Check if the group has a non-zero intensity in WWTP_E
    has.eff = any(treatment == "WWTP_E" & intensity != 0, na.rm = TRUE),

    # Check if the group has a non-zero intensity in the last treatment
    has.last = case_when(
      date == "19.07.2022" & any(treatment == "CMF" & intensity != 0, na.rm = TRUE) ~ TRUE,
      date != "19.07.2022" & any(treatment == "GAC" & intensity != 0, na.rm = TRUE) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  ungroup()

## ============================================== #
## Assign trends - Wide                        ####
## ============================================== #

# Pivot to wide format for trend analysis
df_wide <- df_last_treat %>%
  pivot_wider(names_from = treatment, values_from = intensity, values_fill = 0) %>%
  mutate(
    GAC = ifelse(date == "19.07.2022", NA, GAC),
    O3 = ifelse(date %in% c("15.09.2022", "15.10.2022"), NA, O3),
    AO = ifelse(date %in% c("15.09.2022", "15.10.2022"), AO, NA)
  )

# Flag presence/absence per treatment, treat NA as 0
df_wide_flag <- df_wide %>%
  mutate(
    WWTP_E.0 = coalesce(WWTP_E, 0) != 0,
    O3.0 = coalesce(O3, 0) != 0,
    AO.0 = coalesce(AO, 0) != 0,
    CMF.0 = coalesce(CMF, 0) != 0,
    GAC.0 = coalesce(GAC, 0) != 0,

    # Where is a feature present
    present = purrr::pmap_chr(
      list(WWTP_E.0, O3.0, AO.0, CMF.0, GAC.0),
      ~ {
        idx <- which(c(..1, ..2, ..3, ..4, ..5)) # Which index element is TRUE
        paste(c("WWTP_E", "O3", "AO", "CMF", "GAC")[idx], collapse = "+") # Paste for each TRUE element the treatment
      }
    ),
    present = na_if(present, "")
  ) %>%
  filter(!is.na(present))


# Calculate relative change from WWTP_E to end treatment
df_wide_change <- df_wide_flag %>%
  mutate(
    end_val = coalesce(GAC, CMF, AO, O3, WWTP_E),
    rel_change = 1 - end_val / WWTP_E
  )

# Threshold for trend classification
thr <- 0.35

# Assign trend based on relative change
df_wide_trend <- df_wide_change %>%
  mutate(
    trend = case_when(
      present == first_obs & first_obs != "WWTP_E" &
        is.infinite(rel_change) ~ "NF",
      is.infinite(rel_change) ~ "NF",
      is.na(rel_change) & end_val == 0 ~ "Removed",
      rel_change == 1.0 ~ "Removed",
      abs(rel_change) < thr ~ "Unchanged",
      rel_change > 0 ~ "Decrease",
      rel_change < 0 ~ "Increase",
      TRUE ~ NA_character_
    )
  )

## ============================================== #
## Combine wide dfs                            ####
## ============================================== #

# Combine wide-format data frames for suspects, TPs, and full NTS
df_general_trend_suspects_wide <- df_wide_trend %>%
  filter(group %in% suspect_features$group) %>%
  mutate(set = "Suspect")
df_general_trend_TP_wide <- df_wide_trend %>%
  filter(group %in% TP_features$group) %>%
  mutate(set = "TP")
df_general_trend_suspects_NTS_wide <- df_wide_trend %>% mutate(set = "Full_NTS")
df_general_trend_wide <- rbind(df_general_trend_suspects_wide, df_general_trend_suspects_NTS_wide, df_general_trend_TP_wide) %>% select(-matches(".\\."))


## ============================================== #
## Trends Long - Preparation                   ####
## ============================================== #

# Re-prepare data for long-format trend analysis
df_prepared <- df_nts %>%
  select(group, ret, mz, date, treatment, assigned_pol_mean, intensity)

df_mean_intensity <- df_prepared %>%
  group_by(date, group, assigned_pol_mean, treatment) %>%
  mutate(intensity = mean(intensity, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(date, group, assigned_pol_mean, treatment, .keep_all = TRUE)

df_obs <- df_mean_intensity %>%
  arrange(date, group, assigned_pol_mean, treatment) %>%
  group_by(date, group, assigned_pol_mean) %>%
  mutate(
    present = ifelse(intensity > 0, TRUE, FALSE),
    first_obs = if (any(intensity > 0, na.rm = TRUE)) {
      treatment[intensity > 0][which.min(as.numeric(treatment[intensity > 0]))]
    } else {
      NA
    },
    last_obs = if (any(intensity > 0, na.rm = TRUE)) {
      treatment[intensity > 0][which.max(as.numeric(treatment[intensity > 0]))]
    } else {
      NA
    }
  )


df_last_treat <- df_obs %>%
  group_by(date, group, assigned_pol_mean) %>%
  mutate(
    # Check last treatment step. Depends on the day
    last.treat = last(treatment),

    # Check if the group has a non-zero intensity in WWTP_E
    has.eff = any(treatment == "WWTP_E" & intensity != 0, na.rm = TRUE),

    # Check if the group has a non-zero intensity in the last treatment
    has.last = case_when(
      date == "19.07.2022" & any(treatment == "CMF" & intensity != 0, na.rm = TRUE) ~ TRUE,
      date != "19.07.2022" & any(treatment == "GAC" & intensity != 0, na.rm = TRUE) ~ TRUE,
      TRUE ~ FALSE
    )
  )

df_presence_change <- df_last_treat %>%
  arrange(date, group, assigned_pol_mean, treatment) %>%
  group_by(date, group, assigned_pol_mean) %>%
  mutate(
    # Previous intensity for step-wise change
    prev_intensity = lag(intensity),

    # Relative change per step (e.g., WWTP_E to O3, O3 to CMF, etc.)
    rel_rem_step = 1 - intensity / prev_intensity,

    # Overall relative change (from WWTP_E to current treatment)
    WWTP_E_intensity = first(intensity[treatment == "WWTP_E"]),
    rel_rem_overall = 1 - intensity / WWTP_E_intensity
  )


## ============================================== #
## Assigning Trends - Long                      ####
## ============================================== #

thr <- 0.35

df_general_trend_long <- df_presence_change %>%
  group_by(date, group, assigned_pol_mean) %>%
  mutate(
    # Assign trend
    # Per treatment step
    trend = case_when(
      rel_rem_overall == 1.0 & treatment == "CMF" &
        date == "19.07.2022" & prev_intensity == 0 ~ "Absent",
      rel_rem_step == 1.0 & treatment == "GAC" ~ "Removed",
      rel_rem_step == 1.0 & treatment == "CMF" ~ "Removed",
      rel_rem_overall == 1.0 & rel_rem_step == 1.0 ~ "Removed",
      present == FALSE & is.na(first_obs) & is.na(last_obs) ~ "Absent",
      present == FALSE & is.na(rel_rem_overall) ~ "Absent",
      rel_rem_overall == 1.0 & is.na(rel_rem_step) ~ "Absent",
      first_obs == "WWTP_E" & rel_rem_overall == 0 ~ "start_WWTP_E",
      treatment == first_obs & first_obs != "WWTP_E" &
        is.infinite(rel_rem_overall) ~ "NF",
      abs(rel_rem_step) < thr ~ "Unchanged",
      rel_rem_step > 0 ~ "Decrease",
      is.na(rel_rem_overall) ~ "Unchanged",
      rel_rem_step < 0 ~ "Increase",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

# Combine long-format data frames for suspects, TPs, and full NTS
df_general_trend_suspects_long <- df_general_trend_long %>%
  filter(group %in% suspect_features$group) %>%
  mutate(set = "Suspect")
df_general_trend_TP_long <- df_general_trend_long %>%
  filter(group %in% TP_features$group) %>%
  mutate(set = "TP")
df_general_trend_suspects_NTS_long <- df_general_trend_long %>% mutate(set = "Full_NTS")
df_general_trend_long <- rbind(df_general_trend_suspects_long, df_general_trend_suspects_NTS_long, df_general_trend_TP_long) %>% select(-matches(".\\."))

## =============================================== #
## Specify NF trends further - Long             ####
## =============================================== #

# Further classify newly found (NF) trends
df_NF_trends <- df_general_trend_long %>%
  mutate(trend_s = case_when(
    trend == "Absent" ~ "Absent",
    trend == "start_WWTP_E" ~ "start_WWTP_E",
    trend == "Removed" & is.na(rel_rem_step) ~ "Absent",
    trend %in% c("Increase", "Decrease") &
      is.infinite(rel_rem_step) &
      last_obs == "CMF" &
      date == "19.07.2022" ~ "NF - Pers.",
    trend %in% c("NF") &
      is.infinite(rel_rem_step) &
      last_obs == "CMF" &
      treatment == "CMF" &
      date == "19.07.2022" ~ "NF - Pers.",
    trend %in% c("Increase", "Decrease") &
      is.infinite(rel_rem_step) ~ "NF",
    trend == "Removed" & !is.na(rel_rem_step) ~ "Removed",
    TRUE ~ trend,
    trend == "Removed" ~ paste0("Absent")
  ))

# Identify features with "NF" trends per treatment
df_NF_O3 <- df_general_trend_long %>%
  filter(treatment == "O3" & trend == "NF") %>%
  group_by(group, date) %>%
  select(group, date) %>%
  mutate(group.date = paste0(group, "_", date))
df_NF_CMF <- df_general_trend_long %>%
  filter(treatment == "CMF" & trend == "NF") %>%
  group_by(group, date) %>%
  select(group, date) %>%
  mutate(group.date = paste0(group, "_", date))
df_NF_AO <- df_general_trend_long %>%
  filter(treatment == "AO" & trend == "NF") %>%
  group_by(group, date) %>%
  select(group, date) %>%
  mutate(group.date = paste0(group, "_", date))
df_NF_GAC <- df_general_trend_long %>%
  filter(treatment == "GAC" & trend == "NF") %>%
  group_by(group, date) %>%
  select(group, date) %>%
  mutate(group.date = paste0(group, "_", date))
df_NF_all <- rbind(df_NF_O3, df_NF_CMF, df_NF_AO, df_NF_GAC)

# Refine "NF" trends based on group.date
df_NF_trends <- df_NF_trends %>% mutate(
  group.date = paste0(group, "_", date),
  trend_s = case_when(
    trend_s == "Removed" & treatment == "CMF" & group.date %in% df_NF_all$group.date &
      rel_rem_step == 1.0 ~ "NF - Rem.",
    trend_s == "Removed" & treatment == "GAC" & group.date %in% df_NF_all$group.date &
      rel_rem_step == 1.0 ~ "NF - Rem.",
    trend_s %in% c("Increase", "Decrease", "Unchanged") & treatment == "CMF" & group.date %in% df_NF_all$group.date &
      rel_rem_step != 1.0 ~ "NF - Pers.",
    trend_s %in% c("Increase", "Decrease", "Unchanged") & treatment == "GAC" & group.date %in% df_NF_all$group.date &
      rel_rem_step != 1.0 ~ "NF - Pers.",
    trend_s %in% c("Increase", "Decrease", "Unchanged") ~ "Persistent",
    TRUE ~ trend_s
  )
)

## =============================================== #
## Feature Counts and Intensities               ####
## =============================================== #

# How many features are present per treatment, date, and set
df_feat_n_long <- df_general_trend_long %>%
  filter(intensity != 0) %>%
  group_by(date, treatment, set) %>%
  summarise(n = n(), .groups = "drop")

df_feat_n_wide <- df_general_trend_long %>%
  filter(intensity != 0) %>%
  group_by(date, treatment, set) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = set, values_from = n) %>%
  ungroup()


# What are the total feature intensities per treatment, date, and set #
df_feat_i_tot <- df_general_trend_long %>%
  filter(intensity != 0) %>%
  group_by(date, treatment, set) %>%
  summarise(tot_i = sum(intensity)) %>%
  ungroup()

## =============================================== #
## Trend Analysis - Wide                        ####
## =============================================== #

## How many trends follow each pattern - wide
df_a1 <- df_general_trend_wide %>%
  filter(set == "Suspect") %>%
  group_by(date, trend, set) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(date, set) %>%
  mutate(
    sum_n = sum(n),
    frac = n / sum_n,
    check = sum(frac)
  ) %>%
  arrange(desc(frac)) %>%
  ungroup()

## =============================================== #
## General trends per treatment and set - long  ####
## =============================================== #

# General trends per treatment and set - long
df_b1 <- df_general_trend_long %>%
  filter(treatment != "WWTP_E" & trend != "Absent") %>%
  group_by(treatment, set, trend) %>%
  summarise(n = n(), .groups = "drop") %>%
  ungroup() %>%
  group_by(treatment, set) %>%
  mutate(
    n_total = sum(n),
    frac = n / n_total,
    check = sum(frac)
  )


df_b2 <- df_NF_trends %>%
  filter(treatment != "WWTP_E" & trend_s != "Absent") %>%
  group_by(treatment, set, trend_s) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(treatment, set) %>%
  mutate(
    n_total = sum(n),
    frac = n / n_total,
    check = sum(frac)
  )

## =============================================== #
## Features removed by O3/AO                    ####
## =============================================== #

# Filter for features present in O3/AO and remove last day which contains only data for PFAS CALUX
df_feat_rem_O3AO <- df_general_trend_long %>%
  filter(treatment %in% c("O3", "AO") & trend == "Removed" &
    date != "15.10.2022") %>%
  select(group, date, intensity) %>%
  mutate(presence = ifelse(intensity == 0, 1, 0)) %>%
  pivot_wider(names_from = date, values_from = presence, values_fn = mean) %>%
  select(-intensity) %>%
  mutate(across(where(is.numeric), ~ coalesce(., 0)),
    total_obs = rowSums(across(where(is.numeric), ~ .x == 1))
  )

df_fgroups_rem_O3AO <- df_feat_rem_O3AO %>%
  filter(total_obs == 5) %>%
  pull(group)


## =============================================== #
## Plot trends - wide                           ####
## =============================================== #

# Define trend order and color palettes
trend_order <- rev(c("Removed", "NF", "Unchanged", "Increase", "Decrease"))
trend_order_a2 <- c("NF", "NF - Pers.", "Persistent", "NF - Rem.", "Removed")

shared_colors <- c(
  "Removed" = brewer.pal(n = 8, name = "Paired")[1],
  "Increase" = brewer.pal(n = 8, name = "Paired")[2],
  "NF" = brewer.pal(n = 8, name = "Paired")[3],
  "Decrease" = brewer.pal(n = 8, name = "Paired")[5],
  "Unchanged" = brewer.pal(n = 8, name = "Paired")[6]
)

shared_colors2 <- c(
  "Removed"     = brewer.pal(n = 8, name = "Paired")[1],
  "NF - Rem."   = brewer.pal(n = 8, name = "Paired")[2],
  "NF"          = brewer.pal(n = 8, name = "Paired")[3],
  "NF - Pers."  = brewer.pal(n = 8, name = "Paired")[5],
  "Persistent"  = brewer.pal(n = 8, name = "Paired")[6]
)

# Stacked barplot of trends per date
p1_trends <- df_a1 %>%
  filter(set == "Suspect") %>%
  mutate(
    trend = factor(trend, levels = trend_order),
    frac = frac * 100
  ) %>%
  ggplot(aes(x = date, y = frac, fill = trend, label = n)) +
  geom_bar(stat = "identity", col = "black", position = "stack", linewidth = 0.5) +
  geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = shared_colors) +
  labs(
    x = NULL,
    y = "Relative proportion (%)",
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2))
p1_trends

# Stacked barplot of trends per treatment
set_to_plot <- "Suspect"
p2_trends <- df_b1 %>%
  mutate(
    trend = factor(trend, levels = trend_order),
    frac = frac * 100
  ) %>%
  filter(set == set_to_plot) %>%
  ggplot(
    aes(x = treatment, y = frac, fill = trend, label = n)
  ) +
  geom_bar(stat = "identity", position = "stack", col = "black", linewidth = 0.5) +
  geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = shared_colors) +
  scale_x_discrete(
    labels = c(
      "O3" = expression(paste("O"[3]))
    )
  ) +
  labs(
    x = NULL,
    y = "Relative proportion (%)",
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2))
p2_trends

# Stacked barplot of NF removed and NF persistend trends per treatment
p3_trends <- df_b2 %>%
  mutate(
    trend_s = factor(trend_s, levels = trend_order_a2),
    frac = frac * 100
  ) %>%
  filter(set == set_to_plot) %>%
  ggplot(
    aes(x = treatment, y = frac, fill = trend_s, label = n)
  ) +
  geom_bar(stat = "identity", position = "stack", col = "black", linewidth = 0.5) +
  geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = shared_colors2) + # Use shared colors
  
  scale_x_discrete(
    labels = c(
      "O3" = expression(paste("O"[3]))
    )
  ) +
  labs(
    x = NULL,
    y = "Relative proportion (%)",
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2))
p3_trends


## ============================================== #
# Supplementary Information                    ####
## ============================================== #

## ============================================== #
## Histogram - Removal                          ####
## ============================================== #

p_NTS_histo <- df_general_trend_long %>%
  filter(treatment != "WWTP_E" &
    set %in% c("Full_NTS", "Suspect") &
    trend != "Absent" &
    is.numeric(rel_rem_step)) %>%
  mutate(rel_rem_step = ifelse(is.infinite(rel_rem_step), -1, rel_rem_step)) %>%
  mutate(rel_rem_step = rel_rem_step * 100) %>%
  mutate(set = ifelse(set == "Full_NTS", "All Features", set)) %>%
  ggplot(aes(x = rel_rem_step)) +
  geom_histogram(bins = 100, fill = "lightblue", color = "black") +
  coord_cartesian(xlim = c(-400, NA)) +
  labs(x = "Relative removal (%)", y = "Count") +
  facet_grid(set ~ treatment, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_NTS_histo

## =============================================== #
## Feature removal                               ####
## =============================================== #

p_NTS_boxplot <- df_general_trend_long %>%
  filter(treatment != "WWTP_E" &
    set %in% c("Full_NTS", "Suspect") &
    trend %in% c("Decrease", "Increase", "Unchanged") &
    !is.infinite(rel_rem_step)) %>%
  mutate(rel_rem_step = ifelse(is.infinite(rel_rem_step), -1, rel_rem_step)) %>%
  mutate(rel_rem_step = rel_rem_step * 100) %>%
  mutate(set = ifelse(set == "Full_NTS", "All Features", set)) %>%
  ggplot(aes(x = treatment, y = rel_rem_step, fill = date)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "RdBu") +
  coord_cartesian(ylim = c(-100, 100)) +
  labs(
    x = NULL,
    y = "Relative intensity removal (%)",
    fill = NULL
  ) +
  facet_wrap(. ~ set) +
  theme(legend.position = "bottom")
p_NTS_boxplot

## =============================================== #
## Noise Threshold                              ####
## =============================================== #

# read in raw data ##
raw_data_file_noise_thr <- raw_data_folder[str_detect(raw_data_folder, "noise_thr")]

raw_data_noise_thr <- read.csv(raw_data_file_noise_thr)

# Data wrangling ##
raw_data_noise_thr <- raw_data_noise_thr %>% mutate(par_set = case_when(
  par_set == "a" ~ "Default",
  par_set == "b" & noise_thr == 1000 ~ "Ore. et al. 2025",
  TRUE ~ "This Study"
))

df_long <- raw_data_noise_thr %>%
  pivot_longer(
    cols = c(data_size, total_features, avg_features_p_analysis),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(mode = str_to_title(mode))


# Color palette
my_cols <- brewer.pal(3, "Set1")[c(1, 2)] # red, blue

# individual plots #
p1_total_features <- df_long %>%
  filter(variable == "total_features") %>%
  ggplot(aes(x = noise_thr, y = value, shape = par_set, colour = mode)) +
  geom_point(size = 3) +
  scale_y_continuous(
    labels = label_number_auto(),
    breaks = function(lims) unique(c(pretty_breaks()(lims), 2500000))
  ) +
  scale_shape_discrete(name = "Parameter Set") +
  scale_colour_manual(name = "Ionisation", values = my_cols) +
  labs(
    x = "Noise threshold intensity",
    y = "Total features"
  )
p1_total_features

p2_avg_features <- df_long %>%
  filter(variable == "avg_features_p_analysis") %>%
  ggplot(aes(x = noise_thr, y = value, shape = par_set, colour = mode)) +
  geom_point(size = 3) +
  scale_y_continuous(
    labels = label_number_auto()
  ) +
  scale_shape_discrete(name = "Parameter Set") +
  scale_colour_manual(name = "Ionisation", values = my_cols) +
  labs(
    x = "Noise threshold intensity",
    y = "Features per analysis",
    col = "Ionisation mode"
  )
p2_avg_features

p3_data_size <- df_long %>%
  filter(variable == "data_size") %>%
  ggplot(aes(x = noise_thr, y = value, shape = par_set, colour = mode)) +
  geom_point(size = 3) +
  scale_y_continuous(labels = label_number_auto()) +
  scale_shape_discrete(name = "Parameter set") +
  scale_colour_manual(name = "Ionisation", values = my_cols) +
  labs(
    x = "Noise threshold intensity",
    y = "Data size (MB)",
    col = "Ionisation mode"
  )
p3_data_size

# Add titles
p1_total_features_final <- p1_total_features +
  theme(plot.title = element_text(size = 22, hjust = 0.5, face = "bold"), legend.position = "bottom")

p2_avg_features_final <- p2_avg_features +
  theme(plot.title = element_text(size = 22, hjust = 0.5, face = "bold"), legend.position = "bottom")

p3_data_size_final <- p3_data_size +
  theme(plot.title = element_text(size = 22, hjust = 0.5, face = "bold"), legend.position = "bottom")

## =============================================== #
# Detected feature numbers and intensities      ####
## =============================================== #
p_feat_n <- df_feat_n_long %>%
  mutate(treatment = fct_recode(treatment, "WWTP-E" = "WWTP_E")) %>%
  ggplot(aes(x = treatment, y = n, fill = set)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(
    palette = "Set2",
    labels = c("Full_NTS" = "All Features", "Suspect" = "Suspects", "TP" = "TPs")
  ) +
  labs(
    x = NULL,
    y = expression(n[(Features)]),
    # y = "n(Features)",
    fill = NULL
  ) +
  facet_wrap(. ~ date) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "bottom"
  )
p_feat_n

p_feat_i <- df_feat_i_tot %>%
  mutate(treatment = fct_recode(treatment, "WWTP-E" = "WWTP_E")) %>%
  ggplot(aes(x = treatment, y = tot_i, fill = set)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "",
    y = "Feature intensity"
  ) +
  facet_wrap(. ~ date) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "none"
  )
p_feat_i
