library(dplyr)
library(purrr)
library(readr)
library(lubridate)
library(hms)
library(drc)
library(runner)
library(ggplot2)
library(readxl)

setwd("~/SCIENCE/PhD/5. LittleLeo/Deployments")
# Set your top-level folder here
meta <- read_excel("LittleLeo_MASTER.xlsx")

data_path <- here::here("Data")
processed_path <- here::here("Processed", "01_Processed")
jiggle_path <- here::here("Processed", "02_Jiggle")
speed_path <- here::here("Processed", "03_Speed")
output_path <- speed_path <- here::here("Processed", "04_Merged")


# Compile TBF files #####
# List all matching files recursively
files <- list.files(
  path = data_path,
  pattern = "TBF_Z.csv",
  full.names = TRUE,
  recursive = TRUE
)

# Read, filter, and bind them all
combined_tbf <- files %>%
  map_dfr(~ {
    read_csv(.x, show_col_types = FALSE) %>%
      mutate(source_file = basename(.x),
             ID = substr(source_file, 1, 8)) %>%
      rename(second = RunNo) %>%
      dplyr::select(-source_file)# optional: keep file name
  })


# Compile Depth files
# List all matching files recursively
files <- list.files(
  path = processed_path,
  pattern = "1Hz.csv",
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
    read_csv(.x, col_types = col_spec, show_col_types = FALSE) %>%
      mutate(
        OCDR = as.numeric(OCDR),
        Datetime = as.POSIXct(Datetime.x, format = "%d/%m/%Y %H:%M:%S", tz = "UTC")) %>%
      #Datetime = as.POSIXct(Datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
     dplyr::select(-Datetime.x, -Datetime.y)
  })



# Compile Jiggle files #####
# List all matching files recursively
files <- list.files(
  path = jiggle_path,
  pattern = "jiggle.csv",
  full.names = TRUE,
  recursive = TRUE
)

# Read, filter, and bind them all
combined_jiggle <- files %>%
  map_dfr(~ {
    read_csv(.x, show_col_types = FALSE) %>%
      dplyr::select(Jx_med, Jy_med, Jz_med, Jmag_med, second, ID) # optional: keep file name
  })


# Compile Speed files #####
# List all matching files recursively
files <- list.files(
  path = speed_path,
  pattern = "speed.csv",
  full.names = TRUE,
  recursive = TRUE
)

# Read, filter, and bind them all
combined_speed <- files %>%
  map_dfr(~ {
    read_csv(.x, show_col_types = FALSE) %>%
      dplyr::select(U_pred, U_low, U_high, second, ID) # optional: keep file name
  })

# Compile Video TBF annotation files #####
# List all matching files recursively
files <- list.files(
  path = data_path,
  pattern = "TBF_annotated.csv",
  full.names = TRUE,
  recursive = TRUE
)

# Read, filter, and bind them all
combined_videotbf <- files %>%
  map_dfr(~ {
    read_csv(.x, show_col_types = FALSE) %>%
      filter(Behavior == "Tailbeat") %>%
      filter(is.finite(TBC)) %>%
      mutate(source_file = basename(.x),
             ID = substr(source_file, 1, 8),
             TBF_annotated = 1/TBC) %>%
      dplyr::select(TBF_annotated, TBC, Time, ID) # optional: keep file name
  })

video_start <- combined_d %>%
  filter(Video == 1) %>%
  group_by(ID) %>%
  summarize(
    video_start = first(second)
  ) %>%
  mutate(
    video_start = case_when(video_start < 60*60 ~ 0, # make sure that videos without delay start are not clipped
              TRUE ~ video_start)
  )

combined_videotbf <- combined_videotbf %>%
  left_join(video_start, by = "ID") %>%
  mutate(
    second = floor(Time)+video_start
  ) %>%
  dplyr::select(TBF_annotated, second, ID)

# Calculate mean of annotated TBF as for some reason the annotations have two values per second
combined_videotbf_mean <- combined_videotbf %>%
  group_by(ID, second) %>%
  summarize(
    TBF_annotated = mean(TBF_annotated, na.rm = F)
  ) %>%
  ungroup()

# Check number of observations per deployment
combined_videotbf_mean  %>%
  group_by(ID) %>%
  filter(TBF_annotated > 0.00001) %>%
  summarize(
    n = n()
  ) %>%
  ungroup()


combined_meta_d_jiggle <- combined_d %>%
  left_join(meta %>% 
              rename(ID = DeploymentID) %>%
              dplyr::select(ID, Deployment_start, Sex, TL, Delay_programming)) %>%
  left_join(combined_jiggle, by = c("ID", "second"))

combined_meta_d_jiggle_speed <- combined_meta_d_jiggle %>%
  left_join(combined_speed,
            by = c("ID", "second"))

combined_meta_d_jiggle_speed_tbf <- combined_meta_d_jiggle_speed %>%
  left_join(combined_tbf, by = c("ID", "second"))

combined_meta_d_jiggle_speed_tbf_video <- combined_meta_d_jiggle_speed_tbf %>%
  left_join(combined_videotbf_mean, by = c("ID", "second"), relationship = "one-to-one")



# Check that all numeric columns are filled and none are empty for any deployment
summary_all_numeric <- combined_meta_d_jiggle_speed_tbf_video %>%
  group_by(ID) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )

write.csv(combined_meta_d_jiggle_speed_tbf_video, paste0(output_path, "/All_1Hz_merged_260223.csv"), row.names = F)
