
library(dplyr)
library(zoo)
library(hms)
library(nautilus)
library(data.table)
library(mgcv)

# Define paths ----
data_path <- here::here("Data")
output_path <- here::here("Processed", "01_Processed")

setwd("~/SCIENCE/PhD/5. LittleLeo/Deployments/Data")

# Load calibration file
cal <- read.csv(file.path(data_path, "Calibration/CalFile250806.csv"))

processDive <- function(depth_file, accel_file, start_time, cal) {
  
    
  safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  #----------------------------------------------------------
  # 0. Sanity checks
  #----------------------------------------------------------
  if (!inherits(start_time, "POSIXct")) stop("start_time must be POSIXct")
  if (!file.exists(depth_file)) stop(paste("Depth file not found:", depth_file))
  if (!file.exists(accel_file)) stop(paste("Accel file not found:", accel_file))
  
  message("✅ Starting processing for deployment: ", substr(depth_file, 1, 8))
  
  #----------------------------------------------------------
  # 1. Process depth data
  #----------------------------------------------------------
  message("⏳ Processing depth data...")
  lines_depth <- readLines(depth_file)
  depth <- read.csv(text = lines_depth[6:length(lines_depth)], header = TRUE)
  
  depth <- depth %>%
    mutate(
      RunNo = seq_len(n()),
      Datetime = seq.POSIXt(from = start_time, by = "sec", length.out = n()),
      Time = as_hms(Datetime),
      Depth = Depth - 0.3,
      Depth_3s = rollapply(Depth, width = 3, FUN = mean, align = "center", fill = "extend"),
      Depth_diff = c(NA, diff(Depth_3s))
    )
  
  message("✅ Depth data processed (", nrow(depth), " rows).")
  
  #----------------------------------------------------------
  # 2. Process accelerometer data
  #----------------------------------------------------------
  message("⏳ Processing accelerometer data...")
  lines_acc <- readLines(accel_file)
  start_line <- grep("^ X", lines_acc)
  accel <- read.csv(text = lines_acc[(start_line + 1):length(lines_acc)], header = FALSE)
  colnames(accel) <- c("X", "Y", "Z")
  accel <- accel[, 1:3]
  
  accel <- accel %>%
    mutate(
      ID = substr(depth_file, 1, 8),
      RunNo = seq_len(n()),
      Datetime = seq.POSIXt(from = start_time, by = 0.01, length.out = n()),
      Datetime_str = format(Datetime, "%Y-%m-%d %H:%M:%OS2"),
      Time = as_hms(Datetime),
      Time_str = format(Time, "%H:%M:%OS2")
    )
  
  #----------------------------------------------------------
  # 🧹 Trim start and end where depth < 0.25 continuously
  #----------------------------------------------------------
  message("⏳ Trimming pre- and post-deployment surface intervals...")
  
  # Identify first and last valid depth points
  first_valid <- which(depth$Depth_3s >= 0.25)[1]
  last_valid  <- tail(which(depth$Depth_3s >= 0.25), 1)
  
  if (is.na(first_valid) || is.na(last_valid)) {
    warning("⚠️  No depth >= 0.25 found — skipping trimming.")
  } else {
    # Get time window of valid deployment
    start_valid_time <- depth$Datetime[first_valid]
    end_valid_time   <- depth$Datetime[last_valid]
    
    # Trim depth data
    depth <- depth %>%
      filter(Datetime >= start_valid_time & Datetime <= end_valid_time)
    
    # Trim accelerometer data to same time window
    accel <- accel %>%
      filter(Datetime >= start_valid_time & Datetime <= end_valid_time)
    
    message(paste0("✅ Trimmed data to deployment window: ",
                   format(start_valid_time, "%H:%M:%S"), " — ",
                   format(end_valid_time, "%H:%M:%S")))
  }
  
  
  
  fs <- 100
  accel <- accel %>%
    # remove values > 3 as Little Leonard sensor is only calibrated forvalues below
    mutate(
      across(
        c(X, Y, Z),
        ~ if_else(.x > 3, NA_real_, .x)
      )) %>%
    mutate(
      X = (X - cal$correction_factor[1]) * cal$gain[1],
      Y = (Y - cal$correction_factor[2]) * cal$gain[2],
      Z = (Z - cal$correction_factor[3]) * cal$gain[3]
    ) %>%
    mutate(
      X_stat = rollapply(X, width = 100*3, FUN = mean, align = "center", fill = "extend"),
      Y_stat = rollapply(Y, width = 100*3, FUN = mean, align = "center", fill = "extend"),
      Z_stat = rollapply(Z, width = 100*3, FUN = mean, align = "center", fill = "extend"),
      X_dyna = X - X_stat,
      Y_dyna = Y - Y_stat,
      Z_dyna = Z - Z_stat,
      ODBA = abs(X_dyna) + abs(Y_dyna) + abs(Z_dyna),
      VeDBA = sqrt(X_dyna^2 + Y_dyna^2 + Z_dyna^2),
      jerk_x = c(NA, diff(X)) * fs,
      jerk_y = c(NA, diff(Y)) * fs,
      jerk_z = c(NA, diff(Z)) * fs,
      Jerk = sqrt(jerk_x^2 + jerk_y^2 + jerk_z^2),
      pitch = atan2(Y_stat, sqrt(X_stat^2 + Z_stat^2)),
      roll_raw = atan2(X_stat, Z_stat) # NEW
    )
  
  # save deployment side metadata according to X-axis acceleration
  x_mean <- mean(accel$X_stat, na.rm = TRUE)
  
  deployment_side <- if (x_mean > 0) {
    "left"
  } else if (x_mean < 0) {
    "right"
  } else {
    NA_character_
  }
  
  # correct roll for deployment side
  accel = accel %>%
      mutate(roll = if (x_mean < 0) roll_raw +1.5708 else roll_raw -1.5708)
  
  message("✅ Accelerometer data processed (", nrow(accel), " rows).")
  
  #----------------------------------------------------------
  #  Downsample Accelerometry to 1Hz
  #----------------------------------------------------------
  message("⏳ Downsampling Accelerometry to 1Hz...")
  
  downsample.to = 1
  # convert the desired downsample rate to time interval in seconds
  downsample_interval <- 1 / downsample.to
  
  # select columns to keep
  metrics <- c("X_stat", "Z_stat", "Y_stat", "X_dyna","Z_dyna","Y_dyna", "ODBA", "VeDBA", "Jerk", "roll_raw", "roll", "pitch")
  
  # provide feedback to the user if verbose mode is enabled
  cat("Downsampling data to", downsample.to, "Hz\n")
  
  # round Datetime to the nearest downsample interval
  first_time <- accel$Datetime[1]
  accel1 <- as.data.table(accel)
  accel1[, Datetime := first_time + floor(as.numeric(Datetime - first_time) / downsample_interval) * downsample_interval]
  
  # aggregate numeric metrics using arithmetic mean
  accel1 <- accel1[, lapply(.SD, mean, na.rm=TRUE), by = Datetime, .SDcols = metrics]
  
  accel1 <- accel1 %>%
    mutate(
      ID = accel$ID[1]
    )
  
  # calculate TBF ####
  #tbf <-  calculateTailBeats(accel20, id.col = "ID", datetime.col = "Datetime", motion.col = "Z_dyna", plot.diagnostic = T, plot.output.dir = output_path, n.cores = 3)
  
  #tbf <- tbf[[1]]
  
  message("✅ PDownsampling complete.")
  
  #----------------------------------------------------------
  # 3. Compute pitch offset and OCDR and merge with 1Hz downsampled Acc
  #----------------------------------------------------------
  message("⏳ Calculating pitch offset and OCDR...")
  
  accel1 <- accel1 %>%
    mutate(second = floor(as.numeric(difftime(Datetime, start_time, units = "secs")))) 
  
  depth_accel <- depth %>%
    mutate(second = floor(as.numeric(difftime(Datetime, start_time, units = "secs")))) %>%
    left_join(
      accel1, by = "second"
    ) %>%
    mutate(
      pitch_vv0 = if_else(Depth_diff == 0, pitch, NA_real_),
      pitch_offset_mean = rollapply(pitch_vv0, width = 100*60, FUN = mean, na.rm = TRUE, align = "center", partial = TRUE), # compute pitch correction median with 100 min rolling window
      OCDR = if_else(
        (pitch < -(45 * pi / 180) | pitch > (45 * pi / 180)) & Depth_diff != 0, 
        Depth_diff / sin(pitch),
        NA_real_
      )
    )
  
  # Compute putch correction using GAM model with adaptive smooth and k = 30
  pitch_corr_mod <- gam(pitch_vv0 ~ s(RunNo, bs = "ad", k = 30), method = "REML", data = depth_accel)
  
  # Save model outpur for later
  pitch_corr_fit <- summary(pitch_corr_mod)
  
  # Predict pitch offset using GAM model
  depth_accel$pitch_offset_gam <- predict(pitch_corr_mod, newdata = data.frame(RunNo = depth$RunNo))
  
  depth_accel <- depth_accel %>%
    mutate(pitch_corrected = pitch - pitch_offset_gam)
  
  message("✅ Pitch offset and OCDR calculated.")
  
  #----------------------------------------------------------
  # 4. Map pitch offset & OCDR to 100 Hz timestamps
  #----------------------------------------------------------
  message("⏳ Mapping pitch offset and OCDR to accelerometer data...")
  
  accel <- accel %>%
    mutate(second = floor(as.numeric(difftime(Datetime, start_time, units = "secs")))) %>%
    left_join(
      depth_accel %>%
        select(second, OCDR, Depth, Temp = Temp, pitch_offset_mean, pitch_offset_gam, pitch_corrected),
      by = "second"
    ) 

  
  message("✅ Pitch offset and OCDR mapped to accelerometer data.")
  
  #----------------------------------------------------------
  # 5. Compute summary statistics
  #----------------------------------------------------------
  message("⏳ Computing summary statistics...")
  
  # General depth and temp stats
  depth_stats <- depth_accel %>%
    summarise(
      mean_depth = mean(Depth_3s, na.rm = TRUE),
      max_depth = max(Depth_3s, na.rm = TRUE),
      min_depth = min(Depth_3s, na.rm = TRUE),
      mean_temp = mean(Temp, na.rm = TRUE),
      max_temp = max(Temp, na.rm = TRUE),
      min_temp = min(Temp, na.rm = TRUE),
      mean_pitch_offset = mean(pitch_offset_gam, na.rm = TRUE),
      max_pitch_offset = max(pitch_offset_gam, na.rm = TRUE),
      min_pitch_offset = min(pitch_offset_gam, na.rm = TRUE),
      mean_roll = mean(roll, na.rm = TRUE),
      max_roll = max(roll, na.rm = TRUE),
      min_roll = min(roll, na.rm = TRUE)
      # HERE
    )
  
  # Tagshift
  tagshift <- depth_accel %>%
    summarise(
      start_angle = first(na.omit(pitch_vv0)),
      end_angle   = last(na.omit(pitch_vv0)),
      angle_diff  = {
        s <- first(na.omit(pitch))
        e <- last(na.omit(pitch))
        if (length(s) == 0 || length(e) == 0) NA_real_ else e - s
      }
    )
  
  # Power dives
  depth_accel$power_dive <- abs(depth_accel$pitch) > (45 * pi / 180)
  rle_pd <- rle(depth_accel$power_dive)
  num_powerdives <- sum(rle_pd$values)
  
  powerdives <- depth_accel %>%
    filter(power_dive) %>%
    summarise(
      max_pitch = safe_max(pitch),
      min_pitch = safe_min(pitch),
      max_OCDR = safe_max(abs(OCDR)),
      min_OCDR = safe_min(abs(OCDR)),
      mean_OCDR = mean(abs(OCDR), na.rm = TRUE)
    ) %>%
    mutate(num_dives = num_powerdives)
  
  # Surface
  surface <- depth_accel %>%
    filter(Depth < 2) %>%
    summarise(
      perc_time_surface = n() / nrow(depth) * 100,
      mean_temp_surface = mean(Temp, na.rm = TRUE)
    )
  
  # Build summary row
  
  # to add:
  # min, max, mean, start and end pitch offset
  # min, max, mean roll
  # gam r2
  
  summary_row <- tibble(
    Deployment_ID = substr(depth_file, 1, 8),
    Deployment_start_time = start_valid_time,
    Time_on_animal = difftime(max(depth$Datetime), min(depth$Datetime), units = "hours"),
    Fin_side = deployment_side,
    mean_depth = depth_stats$mean_depth,
    max_depth = depth_stats$max_depth,
    min_depth = depth_stats$min_depth,
    mean_temp = depth_stats$mean_temp,
    max_temp = depth_stats$max_temp,
    min_temp = depth_stats$min_temp,
    pitch_offset_gam_edf = pitch_corr_fit$s.table[1, "edf"],
    pitch_offset_gam_dev_exp = pitch_corr_fit$dev.expl,
    mean_pitch_offset = depth_stats$mean_pitch_offset*180/pi,
    max_pitch_offset = depth_stats$max_pitch_offset*180/pi,
    min_pitch_offset = depth_stats$min_pitch_offset*180/pi,
    mean_roll = depth_stats$mean_roll*180/pi,
    max_roll = depth_stats$max_roll*180/pi,
    min_roll = depth_stats$min_roll*180/pi,
    num_powerdives = powerdives$num_dives,
    max_pitch_powerdive = powerdives$max_pitch*180/pi,
    min_pitch_powerdive = powerdives$min_pitch*180/pi,
    max_OCDR = powerdives$max_OCDR,
    min_OCDR = powerdives$min_OCDR,
    mean_OCDR = powerdives$mean_OCDR,
    start_angle = tagshift$start_angle*180/pi,
    end_angle = tagshift$end_angle*180/pi,
    angle_diff = tagshift$angle_diff*180/pi,
    perc_time_surface = surface$perc_time_surface,
    mean_temp_surface = surface$mean_temp_surface
  )
  
  message("✅ Summary table created.")
  
  #----------------------------------------------------------
  # 6. Final checks
  #----------------------------------------------------------
  if (abs(max(accel$Datetime) - max(depth$Datetime)) > 5) {
    warning("⚠️  Time ranges differ by more than 5 s between files.")
  } else {
    message("✅ Time ranges match within tolerance.")
  }
  
  message("🎉 Processing complete for ", substr(depth_file, 1, 8))
  
  return(list(data1Hz = depth_accel, data100Hz = accel, summary = summary_row, pitch_model = pitch_corr_fit))
}


## LL240719 #####
result <- processDive(
  depth_file = "LL240719/LL240719_D1.txt",
  accel_file = "LL240719/LL240719_A1.txt",
  start_time = as.POSIXct("2024-07-19 14:17:28", tz = "UTC"),
  cal = cal
)

summary <- result$summary

write.csv(result$data1Hz, paste0(output_path, "/LL240719_1Hz.csv"), row.names = F)
write.csv(result$data100Hz, paste0(output_path, "/LL240719_100Hz.csv"), row.names = F)
write.csv(result$summary, paste0(output_path, "/Deployment_summary.csv"), row.names = F)

## LL240825 ####
result <- processDive(
  depth_file = "LL240825/LL240825_D1.txt",
  accel_file = "LL240825/LL240825_A1.txt",
  start_time = as.POSIXct("2024-08-25 22:41:00", tz = "UTC"),
  cal = cal
)

summary <- rbind(summary, result$summary)

write.csv(result$data1Hz, paste0(output_path, "/LL240825_1Hz.csv"), row.names = F)
write.csv(result$data100Hz, paste0(output_path, "/LL240825_100Hz.csv"), row.names = F)
write.csv(result$summary, paste0(output_path, "/Deployment_summary.csv"), row.names = F)

## LL240828 ####
result <- processDive(
  depth_file = "LL240828/LL240828_DT.txt",
  accel_file = "LL240828/LL240828_A.txt",
  start_time = as.POSIXct("2024-08-28 01:27:00", tz = "UTC"),
  cal = cal
)

summary <- rbind(summary, result$summary)

write.csv(result$data1Hz, paste0(output_path, "/LL240828_1Hz.csv"), row.names = F)
write.csv(result$data100Hz, paste0(output_path, "/LL240828_100Hz.csv"), row.names = F)
write.csv(result$summary, paste0(output_path, "/Deployment_summary.csv"), row.names = F)

## LL240920 ####
result <- processDive(
  depth_file = "LL240920/LL240920_DT.txt",
  accel_file = "LL240920/LL240920_A.txt",
  start_time = as.POSIXct("2024-09-21 05:11:33", tz = "UTC"),
  cal = cal
)

summary <- rbind(summary, result$summary)
accel <- result$accel

write.csv(result$data1Hz, paste0(output_path, "/LL240920_1Hz.csv"), row.names = F)
write.csv(result$data100Hz, paste0(output_path, "/LL240920_100Hz.csv"), row.names = F)
write.csv(result$summary, paste0(output_path, "/Deployment_summary.csv"), row.names = F)

## LL250522 ####
result <- processDive(
  depth_file = "LL250522/LL250522_DT.txt",
  accel_file = "LL250522/LL250522_A.txt",
  start_time = as.POSIXct("2025-05-23 12:22:28", tz = "UTC"),
  cal = cal
)

summary <- rbind(summary, result$summary)
accel <- result$accel

write.csv(result$data1Hz, paste0(output_path, "/LL250522_1Hz.csv"), row.names = F)
write.csv(result$data100Hz, paste0(output_path, "/LL250522_100Hz.csv"), row.names = F)
write.csv(result$summary, paste0(output_path, "/Deployment_summary.csv"), row.names = F)

## LL250709 ####
result <- processDive(
  depth_file = "LL250709/LL250709_DT.txt",
  accel_file = "LL250709/LL250709_A.txt",
  start_time = as.POSIXct("2025-07-09 23:50:02", tz = "UTC"),
  cal = cal
)

summary <- rbind(summary, result$summary)
accel <- result$accel

write.csv(result$data1Hz, paste0(output_path, "/LL250709_1Hz.csv"), row.names = F)
write.csv(result$data100Hz, paste0(output_path, "/LL250709_100Hz.csv"), row.names = F)
write.csv(result$summary, paste0(output_path, "/Deployment_summary.csv"), row.names = F)

## LL250721 ####
result <- processDive(
  depth_file = "LL250721/LL250721_DT.txt",
  accel_file = "LL250721/LL250721_A.txt",
  start_time = as.POSIXct("2025-07-21 22:33:00", tz = "UTC"),
  cal = cal
)

summary <- rbind(summary, result$summary)

write.csv(result$data1Hz, paste0(output_path, "/LL250721_1Hz.csv"), row.names = F)
write.csv(result$data100Hz, paste0(output_path, "/LL250721_100Hz.csv"), row.names = F)
write.csv(result$summary, paste0(output_path, "/Deployment_summary.csv"), row.names = F)

## LL250915 ####
result <- processDive(
  depth_file = "LL250915/LL250915_DT.txt",
  accel_file = "LL250915/LL250915_A.txt",
  start_time = as.POSIXct("2025-09-16 23:19:58", tz = "UTC"),
  cal = cal
)

summary <- rbind(summary, result$summary)

write.csv(result$data1Hz, paste0(output_path, "/LL250915_1Hz.csv"), row.names = F)
write.csv(result$data100Hz, paste0(output_path, "/LL250915_100Hz.csv"), row.names = F)
write.csv(summary, paste0(output_path, "/Deployment_summary.csv"), row.names = F)
