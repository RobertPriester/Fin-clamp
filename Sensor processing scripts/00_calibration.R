# ============================================================
#  Project: LL
#  Script:  Calibration
#  Purpose: Processing calibration files and calculating calibration metrics
#  Author:  Robert Priester 
#  Date:    251106
#  Version: [v01]
#
#  R. File naming system: [project_acronym]_[step]_[content]_[yymmdd].R
#
# ============================================================
#  Description:
#  [Brief paragraph describing purpose, key inputs, outputs, 
#   and expected workflow position (e.g., "Step 2: data cleaning 
#   before analysis").]
#
#  Input files:
#  - 250806_CAL_procc.csv  : processed accelerometry file
#
#  Output files:
#  - [output1.csv]    : [description]
#  - [output2.png]    : [description]
#
#  Dependencies:
#  - R version: [e.g. 4.4.1]
#  - Required packages: [dplyr, ggplot2, mgcv, etc.]
#
#  Notes:
#  [Special remarks, assumptions, or limitations]
# ============================================================

# Load packages ----
library(dplyr)
library(ggplot2)
library(tidyr)
library(here)

# Define paths ----
data_path <- here::here("Data/", "Calibration", "250806_CAL")
output_path <- here::here("Data/", "Calibration", "250806_CAL")

# Load data ----
data <- read.csv(file.path(data_path, "250806_CAL_procc.csv"))

steady <- data[501: nrow(data),] %>%
  filter(ODBA < 0.03, Jerk < 10)

baseline_table <- steady %>%
  pivot_longer(cols = c(X, Y, Z), names_to = "axis", values_to = "value") %>%
  filter(between(value, -0.1, 0.1)) %>%
  group_by(axis) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    .groups = "drop"
  )

steady_baseline_corr <- steady %>%
  mutate(
    X = X - baseline_table$median[1],
    Y = Y - baseline_table$median[2],
    Z = Z - baseline_table$median[3]
  )


correction_positive <- steady %>%
  pivot_longer(cols = c(X, Y, Z), names_to = "axis", values_to = "value") %>%
  filter(between(value, 0.9, 1.1)) %>%
  group_by(axis) %>%
  summarise(
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

correction_negative <- steady %>%
  pivot_longer(cols = c(X, Y, Z), names_to = "axis", values_to = "value") %>%
  filter(between(value, -1.1, -0.9)) %>%
  group_by(axis) %>%
  summarise(
    min = min(value, na.rm = TRUE),
    .groups = "drop"
  )


calibration_file <- cbind(baseline_table, correction_negative[,2], correction_positive[,2]) %>%
  rename(
    offset = median
  ) %>%
  mutate(
    correction_factor = (min + max)/2,
    gain = 1/(max-correction_factor)
  )



# testing calibration
calibrated_file <- steady%>%
  mutate(
    X = X-calibration_file$correction_factor[1],
    Y = Y-calibration_file$correction_factor[2],
    Z = Z-calibration_file$correction_factor[3]
    ) %>%
  mutate(
    X = X*calibration_file$gain[1],
    Y = Y*calibration_file$gain[2],
    Z = Z*calibration_file$gain[3]
  )

max(calibrated_file$Y)

# Analysis ----
# [your code here]

# Save outputs ----
write.csv(calibration_file, file.path(output_path, "CalFile250806.csv"), row.names = F)

# Session info ----
sessionInfo()