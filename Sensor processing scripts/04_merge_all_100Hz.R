library(dplyr)
library(purrr)
library(readr)
library(lubridate)
library(hms)
library(drc)
library(runner)
library(ggplot2)
library(readxl)
library(fst)
library(data.table)
library(tidyr)

setwd("~/SCIENCE/PhD/5. LittleLeo/Deployments")
# Set your top-level folder here
meta <- read_excel("LittleLeo_MASTER.xlsx")

data_path <- here::here("Data")
processed_path <- here::here("Processed", "01_Processed")

output_path <- here::here("Processed", "04_Merged")


# Compile Depth files
# List all matching files recursively
files <- list.files(
  path = processed_path,
  pattern = "100Hz.csv",
  full.names = TRUE,
  recursive = TRUE
)

col_spec <- cols(
  OCDR = col_double(),
  .default = col_guess()
)

# Read, filter, and bind them all
combined_d <- files %>%
  map_dfr(~ {
    read_csv(.x, col_types = col_spec, show_col_types = FALSE) #%>%
      #mutate(
        #OCDR = as.numeric(OCDR),
        #Datetime = as.POSIXct(Datetime.x, format = "%d/%m/%Y %H:%M:%S", tz = "UTC")) %>%
      #Datetime = as.POSIXct(Datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
      #dplyr::select(-Datetime.x, -Datetime.y)
  })

combined_meta_d <- combined_d %>%
  left_join(meta %>% 
              rename(ID = DeploymentID) %>%
              dplyr::select(ID, Deployment_start, Sex, TL, Delay_programming))

write_fst(combined_meta_d, paste0(output_path, "/All_100Hz_merged_260206.fst"), compress = 100)

setDT(combined_meta_d)

summary <- combined_meta_d[
  ,
  lapply(.SD, function(x)
    list(
      min  = min(x, na.rm = TRUE),
      mean = mean(x, na.rm = TRUE),
      max  = max(x, na.rm = TRUE)
    )
  ),
  by = ID,
  .SDcols = c("X_dyna", "Y_dyna", "Z_dyna", "VeDBA", "jerk_x", "jerk_y", "jerk_z", "Jerk")
]

summary <- summary[
  ,
  lapply(.SD, unlist),
  by = ID
]

summary_clean <- summary %>%
  filter(
    is.finite(X_dyna)
  )%>%
  mutate(statistic = rep(c("min", "mean", "max"), length.out = n()))


summary_wide <- summary_clean %>%
  pivot_wider(
    names_from  = statistic,
    values_from = where(is.numeric),
    names_glue  = "{.value}_{statistic}"
  )


write.csv(summary_wide, paste0(output_path, "/All_100Hz_summary.csv"), row.names = F)
