library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(readxl)


## ============================================= #
# Read in data                                ####
## ============================================= #

# List all files in the  directory, read in data
raw_data_folder <- list.files("raw_data/EBM", full.names = TRUE)
raw_data_file_EBM <- raw_data_folder[str_detect(raw_data_folder, "raw_data_EBM.xlsx")]
raw_data_EBM <- read_excel(raw_data_file_EBM)

# assign to working df, convert to lowercase
df_raw_EBM <- raw_data_EBM %>% rename_with(str_to_lower)

## ============================================= #
# Data Wrangling                              ####
## ============================================= #

# Lookup tables for dates
lookup_dates <- data.frame(
  day = c(paste0("D", 1:6)),
  date = paste0(c("19.07.", "26.07.", "16.08.", "26.08.", "15.09.", "15.10."), "2022")
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
    below_loq = ifelse(str_detect(result, "LOQ|< "), 1, 0),
    result = as.numeric(result),
    result = ifelse(is.na(result), 0, result),
    result = ifelse(unit == "ug/l", result * 1000, result),
    unit = ifelse(unit == "ug/l", "ng/l", unit),
    result_loq = as.numeric(ifelse(below_loq == 1, loq, result))
  ) %>%
  filter(date %in% lookup_dates$date)

# Define order and change sample/treatment names
treat_order <- c("WWTP-E", "O3", "AO", "CMF", "GAC")
df_EBM <- df_EBM %>%
  rename(treatment = sample) %>%
  mutate(
    treatment = factor(
      treatment,
      levels = treat_order
    ),
    date = factor(date, levels = c(lookup_dates$date))
  )

# Generate df for PFAS-CALUX
df_PFAS <- df_EBM %>% filter(calux_assay == "PFAS")

# Remove outlier PFAS-CALUX
outlier_pos <- df_EBM %>%
  filter(calux_assay == "PFAS") %>%
  slice_max(result_loq, n = 1, with_ties = FALSE)
df_EBM <- df_EBM %>% filter(rowid != outlier_pos$rowid)

# Calculate mean effect per bioassay, day, and treatment
df_EBM_mean_eff <- df_EBM %>%
  group_by(treatment, calux_assay, date) %>%
  summarise(
    mean_eff = mean(result, na.rm = TRUE),
    sd_eff = sd(result, na.rm = TRUE),
    mean_loq_eff = mean(result_loq, na.rm = TRUE),
    sd_loq_eff = sd(mean_loq_eff, na.rm = TRUE)
  ) %>%
  mutate(unit = "ng/L") %>%
  left_join(lookup_ebt, join_by(calux_assay))

# Generate a list of all assays for further processing
df_list <- df_EBM_mean_eff %>%
  filter(!is.na(treatment)) %>% # filter NAs due to FB, Lab Blanks etc.
  group_by(calux_assay) %>%
  group_split()


# Calculate mean removal per list element
df_list_output <- df_list %>%
  map(~ .x %>%
    group_by(date) %>%
    mutate(
      baseline            = mean_loq_eff[treatment == "WWTP-E"],
      rel_removal_overall = (baseline - mean_loq_eff) / baseline * 100,
      rel_removal_overall = ifelse(rel_removal_overall == 0, NA, rel_removal_overall),
      rel_removal_overall = round(rel_removal_overall)
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
df_EBM_mean_eff <- df_EBM_mean_eff %>% mutate(
  
  # Surface water exceedance
  ebt_sw_ex = mean_loq_eff / ebt_sw_ngl,
  ebt_sw_df10_ex = mean_loq_eff / ebt_sw_df10,
  
  # Drinking water exceedance
  ebt_dw_ex = mean_loq_eff / ebt_dw_ngl,
  ebt_dw_df10_ex = mean_loq_eff / ebt_dw_df10,
  
  # Categorise exceedance for surface water
  ebt_sw_ex_cat = case_when(
    ebt_sw_ex <= 1 ~ "C0",
    ebt_sw_ex > 1 & ebt_sw_ex <= 3 ~ "C1",
    ebt_sw_ex > 3 & ebt_sw_ex <= 10 ~ "C2",
    ebt_sw_ex > 10 & ebt_sw_ex <= 100 ~ "C3",
    ebt_sw_ex > 100 ~ "C4",
    TRUE ~ NA_character_
  ),
  ebt_sw_df10_ex_cat = case_when(
    ebt_sw_df10_ex <= 1 ~ "C0",
    ebt_sw_df10_ex > 1 & ebt_sw_df10_ex <= 3 ~ "C1",
    ebt_sw_df10_ex > 3 & ebt_sw_df10_ex <= 10 ~ "C2",
    ebt_sw_df10_ex > 10 & ebt_sw_df10_ex <= 100 ~ "C3",
    ebt_sw_df10_ex > 100 ~ "C4",
    TRUE ~ NA_character_
  ),
  
  # Categorise exceedance for drinking water
  ebt_dw_ex_cat = case_when(
    ebt_dw_ex <= 1 ~ "C0",
    ebt_dw_ex > 1 & ebt_dw_ex <= 3 ~ "C1",
    ebt_dw_ex > 3 & ebt_dw_ex <= 10 ~ "C2",
    ebt_dw_ex > 10 & ebt_dw_ex <= 100 ~ "C3",
    ebt_dw_ex > 100 ~ "C4",
    TRUE ~ NA_character_
  ),
  ebt_dw_df10_ex_cat = case_when(
    ebt_dw_df10_ex <= 1 ~ "C0",
    ebt_dw_df10_ex > 1 & ebt_dw_df10_ex <= 3 ~ "C1",
    ebt_dw_df10_ex > 3 & ebt_dw_df10_ex <= 10 ~ "C2",
    ebt_dw_df10_ex > 10 & ebt_dw_df10_ex <= 100 ~ "C3",
    ebt_dw_df10_ex > 100 ~ "C4",
    TRUE ~ NA_character_
  )
)


## ============================================= #
# Plots                                       ####
## ============================================= #

## ============================================= #
## Detected Effects                           ####
## ============================================= #



# Change level (plot) order
df_EBM_mean_eff <- df_EBM_mean_eff %>%
  mutate(calux_assay = factor(calux_assay, levels = c(
    "PXR", "Nrf2", "Cytotox", "PAH", "ERα", "PFAS"
  )))

## ============================================= #
### BEQ (Biological Equivalent Concentration) ####
## ============================================= #

pos <- position_dodge(width = 0.3)
p1_BEQ <- df_EBM_mean_eff %>%
  filter(!is.na(treatment)) %>%
  mutate(mean_loq_eff = mean_loq_eff / 1000) %>% # transform to µg/L
  ggplot(aes(x = treatment, y = mean_loq_eff, fill = date, group = date, col = date)) +
  geom_point(
    size = 5,
    pch = 21,
    position = pos,
    na.rm = TRUE,
    alpha = 0.8,
    col = "black"
  ) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(
    "O3" = expression(paste("O"[3]))
  )) +
  labs(
    x = NULL,
    y = expression(paste("BEQ (µgL"^-1, ")")),
    fill = NULL
  ) +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold"),
    axis.text.x = element_text(angle = 40, hjust = 1)
  ) +
  facet_wrap(. ~ calux_assay, scales = "free")
p1_BEQ

## ============================================= #
### Effect removal                            ####
## ============================================= #

# Generate every plot individually to control axis
my_cols <- brewer.pal(6, "Set2")
assay_list <- levels(df_EBM_mean_eff$calux_assay)

y_ranges <- list(
  "PXR"     = c(0, 100),
  "PAH"     = c(50, 100),
  "ERα"     = c(50, 100),
  "Cytotox" = c(-160, 100),
  "Nrf2"    = c(-160, 100),
  "PFAS"    = c(0, 100)
)
y_breaks <- list(
  "PXR"     = seq(0, 100, by = 25),
  "PAH"     = seq(50, 100, by = 25),
  "ERα"     = seq(50, 100, by = 25),
  "Cytotox" = seq(-150, 100, by = 50),
  "Nrf2"    = seq(-150, 100, by = 50),
  "PFAS"    = seq(0, 100, by = 25)
)

# Define date colors for plotting
date_levels <- levels(df_EBM_mean_eff$date)
date_cols <- brewer.pal(max(3, length(date_levels)), "Set2")[seq_along(date_levels)]
names(date_cols) <- date_levels

# Generate plots for each assay
all_plots <- imap(assay_list, ~ {
  assay <- .x

  df_plot <- df_EBM_mean_eff %>%
    filter(treatment != "WWTP-E", calux_assay == assay)

  ggplot(df_plot, aes(x = treatment, y = rel_removal_overall, group = date)) +
    geom_point(
      aes(fill = date),
      position = pos,
      size = 5,
      shape = 21,
      colour = "black",
      show.legend = TRUE
    ) +
    scale_fill_manual(values = date_cols, drop = FALSE) +
    scale_y_continuous(limits = y_ranges[[assay]], breaks = y_breaks[[assay]]) +
    scale_x_discrete(labels = c(
      "O3" = expression(paste("O"[3]))
    )) +
    labs(x = "", y = "Effect removal (%)", fill = "") +
    facet_wrap(~calux_assay) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "white"),
      strip.text.x = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
})

# Combine plots into a single layout
p2_BEQ_Eff <- wrap_plots(all_plots, ncol = 3) +
  plot_layout(axis_titles = "collect", guides = "collect")
p2_BEQ_Eff


## ============================================= #
# Risk Analysis                               ####
## ============================================= #

## ============================================= #
## EBTs (Effect-Based Triggers)               ####
## ============================================= #

# Define alpha value for transparency in background bands
# and dodging position
alpha_value <- 0.2
pos2 <- position_dodge2(width = 0.4)

# Define treatment levels and colors
treat_levels <- levels(df_EBM_mean_eff$treatment)
treat_cols <- brewer.pal(max(3, length(treat_levels)), "Set2")[seq_along(treat_levels)]
names(treat_cols) <- treat_levels

# Define bands for EBT exceedance categories
band_labels <- tibble::tribble(
  ~label, ~ymin, ~ymax,
  "C[0]", 0.45, 1,
  "C[1]", 1, 3,
  "C[2]", 3, 10,
  "C[3]", 10, 100,
  "C[4]", 100, 200
)

# Calculate midpoint for band labels
band_labels <- band_labels %>%
  mutate(mid = sqrt(ymin * ymax))

# Plot EBT exceedance for surface water
p1_point_EBT_SW <- df_EBM_mean_eff %>%
  filter(!is.na(ebt_sw_ngl)) %>%
  ggplot(aes(x = calux_assay, y = ebt_sw_ex)) +
  annotate("rect",
    fill = "white", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 0, ymax = 1
  ) +
  annotate("rect",
    fill = "green", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 1, ymax = 3
  ) +
  annotate("rect",
    fill = "yellow", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 3, ymax = 10
  ) +
  annotate("rect",
    fill = "orange", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 10, ymax = 100
  ) +
  geom_point(aes(fill = treatment),
    position = pos2,
    size = 6,
    shape = 21,
    colour = "black"
  ) +
  geom_text(
    data = filter(band_labels, label != "C[4]"),
    aes(x = Inf, y = mid, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    fontface = "bold",
    size = 8,
    parse = TRUE # To parse labels with subscripts
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
  scale_y_log10() +
  theme(legend.position = "bottom")
p1_point_EBT_SW

# Plot EBT exceedance for drinking water
p2_point_EBT_DW <- df_EBM_mean_eff %>%
  filter(!is.na(ebt_dw_ngl)) %>%
  ggplot(aes(x = calux_assay, y = ebt_dw_ex)) +
  annotate("rect",
    fill = "white", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 1, ymax = 3
  ) +
  annotate("rect",
    fill = "green", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 1, ymax = 3
  ) +
  annotate("rect",
    fill = "yellow", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 3, ymax = 10
  ) +
  annotate("rect",
    fill = "orange", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 10, ymax = 100
  ) +
  annotate("rect",
    fill = "red", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 100, ymax = Inf
  ) +
  geom_point(aes(fill = treatment),
    position = pos2,
    size = 6,
    shape = 21,
    colour = "black"
  ) +
  scale_y_log10() +
  geom_text(
    data = band_labels,
    aes(x = Inf, y = mid, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    fontface = "bold",
    size = 8,
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
  theme(
    legend.position = "none"
  ) +
  scale_fill_manual(values = treat_cols, drop = TRUE)
p2_point_EBT_DW

## ============================================= #
## EBTs - With DF10                           ####
## ============================================= #

# Prepare data for surface water EBT exceedance with DF10
df_EBT_SW <- df_EBM_mean_eff %>%
  filter(!is.na(ebt_sw_ngl)) %>%
  select(treatment, calux_assay, ebt_sw_ex)
df_EBT_SW_DF10_only <- df_EBM_mean_eff %>%
  filter(!is.na(ebt_sw_ngl)) %>%
  select(treatment, calux_assay, ebt_sw_df10_ex) %>%
  mutate(calux_assay = paste0(calux_assay, "_DF10"))
df_EBT_SW_DF10 <- rbind(df_EBT_SW, df_EBT_SW_DF10_only)

# Plot for surface water EBT exceedance with DF10
pDF10SW <- df_EBT_SW_DF10 %>%
  ggplot(aes(x = calux_assay)) +
  annotate("rect",
    fill = "white", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 0, ymax = 1
  ) +
  annotate("rect",
    fill = "green", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 1, ymax = 3
  ) +
  annotate("rect",
    fill = "yellow", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 3, ymax = 10
  ) +
  annotate("rect",
    fill = "orange", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 10, ymax = 100
  ) +
  geom_point(
    aes(
      x = calux_assay,
      y = ebt_sw_ex,
      fill = treatment
    ),
    position = pos2,
    size = 5,
    shape = 21,
    colour = "black"
  ) +
  geom_point(
    aes(
      x = calux_assay,
      y = ebt_sw_df10_ex,
      fill = treatment
    ),
    position = pos2,
    size = 5,
    shape = 21,
    colour = "black"
  ) +
  geom_text(
    data = filter(band_labels, label != "C[4]"),
    aes(x = Inf, y = mid, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    fontface = "bold",
    size = 4,
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
  scale_x_discrete(
    labels = c(
      "ERα_DF10" = expression(paste("ERα"[DF10])),
      "Nrf2_DF10" = expression(paste("Nrf2"[DF10])),
      "PXR_DF10" = expression(paste("PXR"[DF10])),
      "PAH_DF10" = expression(paste("PAH"[DF10]))
    )
  ) +
  scale_y_log10() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "bottom" 
  )
pDF10SW


# Prepare data for drinking water EBT exceedance with DF10
df_EBT_DW <- df_EBM_mean_eff %>%
  filter(!is.na(ebt_dw_ngl)) %>%
  select(treatment, calux_assay, ebt_dw_ex)
df_EBT_DW_DF10_only <- df_EBM_mean_eff %>%
  filter(!is.na(ebt_dw_ngl)) %>%
  select(treatment, calux_assay, ebt_dw_df10_ex) %>%
  mutate(calux_assay = paste0(calux_assay, "_DF10"))
df_EBT_DW_DF10 <- rbind(df_EBT_DW, df_EBT_DW_DF10_only)

# Plot for drinking water EBT exceedance with DF10
pDF10DW <- df_EBT_DW_DF10 %>%
  ggplot(aes(x = calux_assay)) +
  annotate("rect",
    fill = "white", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 0, ymax = 1
  ) +
  annotate("rect",
    fill = "green", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 1, ymax = 3
  ) +
  annotate("rect",
    fill = "yellow", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 3, ymax = 10
  ) +
  annotate("rect",
    fill = "orange", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 10, ymax = 100
  ) +
  annotate("rect",
    fill = "red", alpha = alpha_value,
    xmin = -Inf, xmax = Inf,
    ymin = 100, ymax = Inf
  ) +
  geom_point(
    aes(
      x = calux_assay,
      y = ebt_dw_ex,
      fill = treatment
    ),
    position = pos2,
    size = 5,
    shape = 21,
    colour = "black"
  ) +
  geom_point(
    aes(
      x = calux_assay,
      y = ebt_dw_df10_ex,
      fill = treatment
    ),
    position = pos2,
    size = 5,
    shape = 21,
    colour = "black"
  ) +
  geom_text(
    data = band_labels,
    aes(x = Inf, y = mid, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    fontface = "bold",
    size = 4,
    parse = TRUE
  ) +
  coord_cartesian(
    clip = "off" # drawing outside panel
  ) +
  labs(
    x = "",
    y = "EBT-Exceedance",
    fill = ""
  ) +
  scale_fill_manual(values = treat_cols, drop = TRUE) +
  scale_x_discrete(
    labels = c(
      "ERα_DF10" = expression(paste("ERα"[DF10])),
      "Nrf2_DF10" = expression(paste("Nrf2"[DF10])),
      "PFAS_DF10" = expression(paste("PFAS"[DF10])),
      "PAH_DF10" = expression(paste("PAH"[DF10]))
    )
  ) +
  scale_y_log10() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "none"
  )
pDF10DW