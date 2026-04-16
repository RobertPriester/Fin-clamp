library(dplyr)
library(zoo)
library(hms)
library(dplyr)
library(readr)
library(purrr)


setwd("~/SCIENCE/PhD/5. LittleLeo/Deployments")
# Set your top-level folder here

data_path <- here::here("Processed", "01_Merged")
folder_path = data_path

output_path <- here::here("Processed", "02_Jiggle")
output_path2 <- here::here("Processed", "02_Jiggle", "Model_fit")


# List all matching files recursively
files <- list.files(
  path = folder_path,
  pattern = "_100Hz\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

o <- read.csv(files[8]) # cycle manually through all files


library(data.table); library(signal); library(zoo); library(ggplot2)

# === inputs & parameters ====================================================
# assume your dataframe is named `a` (as shown in head(a))
DT <- as.data.table(o)

fs <- 100                     # accelerometer sampling rate (Hz)
bp_low <- 10                  # lower bandpass cutoff (Hz)
bp_high <- 45                 # upper bandpass cutoff (Hz) -- must be < fs/2
win_sec <- 0.5                # RMS window length (seconds) used for initial RMS (centered)
agg_bin <- 1                  # aggregate bin for regression (seconds) - usually 1s (OCDR is 1Hz)
pitch_thresh_rad <- pi/4      # 45 degrees in radians (for selecting power dives)
min_OCDR_for_fit <- 0.01      # remove zeros/near-zero speeds prior to log regression


# === Sanity checks ==========================================================
if(bp_high >= fs/2) stop("bp_high must be less than Nyquist (fs/2). Choose a smaller bp_high.")
if(!all(c("X","Y","Z","Time","second","pitch_corrected","OCDR") %in% names(DT))){
  stop("Required columns missing. Need at least X,Y,Z,Time,second,pitch_corrected,OCDR.")
}

# convert Time to POSIXct if needed (optional)
# DT[, Time := as.POSIXct(Time, format="%H:%M:%OS", tz="UTC")]  # adapt to your time format if needed

# === 1) design bandpass FIR filter (symmetric) =================================
# We'll use a linear-phase FIR filter and filtfilt to avoid phase shift.
# Choose filter order: use 127 (128-point symmetric FIR ~ order 127)
filt_order <- 127
Wn <- c(bp_low, bp_high) / (fs/2)   # normalized cutoffs
bp_filt <- fir1(n = filt_order, w = Wn, type = "pass", window = hamming(filt_order+1))

# helper: zero-phase filtering (filtfilt)
filtfilt_safe <- function(b, x){
  # use signal::filtfilt; if numerical issues, fallback to forward/backward filtfilt
  as.numeric(filtfilt(b, a = 1, x))
}

# === 2) bandpass each axis and form high-frequency magnitude ==================
# Use columns X,Y,Z as the raw accel (units in g) — adjust names if different
DT[, X_bp := filtfilt_safe(bp_filt, X)]
DT[, Y_bp := filtfilt_safe(bp_filt, Y)]
DT[, Z_bp := filtfilt_safe(bp_filt, Z)]

# magnitude of bandpassed (high-frequency) vector
DT[, mag_bp := sqrt(X_bp^2 + Y_bp^2 + Z_bp^2)]

# === 3) compute RMS in centered sliding window (win_sec) =======================
win_n <- round(win_sec * fs)        # samples per window (e.g., 0.5s * 100 = 50 samples)
if(win_n %% 2 == 0) win_n <- win_n + 1  # make odd for a centered window (optional)

# RMS = sqrt(mean(mag_bp^2)) over window
DT[, jig_rms := sqrt( rollapply(mag_bp^2, width = win_n, FUN = mean, fill = NA, align = "center") )]

# === 4) convert RMS to dB (amplitude dB) =====================================
# Use 20*log10 for amplitude-to-dB. Add small epsilon to avoid -Inf.
eps <- .Machine$double.eps * 1e3
DT[, jig_dB := 20 * log10(pmax(jig_rms, eps))]

# === 5) aggregate to 1 Hz bins (match OCDR sampling) ==========================
# We assume `second` column in DT indexes seconds (integer); if not, create it
# If your `second` column is already present and matches time, use it.
if(!"second" %in% names(DT)){
  # create second index relative to first timestamp
  DT[, second := floor(as.numeric(difftime(Time, Time[1], units="secs")))]
}

# compute per-second summaries (median is robust)
per_sec <- DT[, .(
  time = first(Time),
  Jx_med = median(20*log10(pmax(sqrt( rollapply(X_bp^2, width = win_n, FUN = mean, fill = NA, align = "center") ), eps)), na.rm=TRUE),
  Jy_med = median(20*log10(pmax(sqrt( rollapply(Y_bp^2, width = win_n, FUN = mean, fill = NA, align = "center") ), eps)), na.rm=TRUE),
  Jz_med = median(20*log10(pmax(sqrt( rollapply(Z_bp^2, width = win_n, FUN = mean, fill = NA, align = "center") ), eps)), na.rm=TRUE),
  Jmag_med = median(jig_dB, na.rm=TRUE),
  pitch_med = median(pitch_corrected, na.rm=TRUE),
  OCDR_med = median(OCDR, na.rm=TRUE),
  Depth_med = median(Depth, na.rm=TRUE)
), by = second]

# === 6) select power-dive seconds (|pitch| > 45 deg) ===========================
per_sec[, pitch_deg := pitch_med * 180 / pi]
power_secs <- per_sec[abs(pitch_med) > pitch_thresh_rad & !is.na(OCDR_med)]

cat("Number of 1-s bins with |pitch|>45deg:", nrow(power_secs), "\n")
# You reported 39 power dives — make sure this equals nrow(power_secs) or adjust by event grouping.

# === 7) prepare OCDR and jiggle for regression =================================
# OCDR must be positive - use absolute value if needed (descending vs ascending signs)
power_secs <- power_secs[!is.na(OCDR_med) & (abs(OCDR_med) > min_OCDR_for_fit)]

# optionally take absolute OCDR or only positive OCDR depending on sign convention
power_secs[, OCDR_abs := abs(OCDR_med)]

# check distributions
summary(power_secs[, .(Jx_med, Jy_med, Jz_med, Jmag_med, OCDR_abs)])
hist(power_secs$Jmag_med, main="J (dB) for power-dive seconds", xlab="J (dB)")

# === 8) fit exponential relationship: log(OCDR) ~ linear(Jx,Jy,Jz) ============
# This is equivalent to OCDR = a * exp(b1*Jx + b2*Jy + b3*Jz)
fit_df <- power_secs[OCDR_abs > min_OCDR_for_fit]
fit_df[, logO := log(OCDR_abs)]

clean_fit <- fit_df %>%
  dplyr::filter(Jmag_med < -10 & OCDR_abs > 0.5)

# Fit linear model on log scale (quick and robust)
lm_fit <- lm(logO ~ Jx_med + Jy_med + Jz_med, data = clean_fit)
summary(lm_fit)

# Convert coefficients back to exponential form for prediction:
# logO = beta0 + beta1*Jx + beta2*Jy + beta3*Jz
# => OCDR_pred = exp(beta0 + sum(beta_i * Ji))
fit_df[, OCDR_pred := exp(predict(lm_fit, newdata = fit_df))]

# goodness of fit
rss <- sum((fit_df$OCDR_abs - fit_df$OCDR_pred)^2)
tss <- sum((fit_df$OCDR_abs - mean(fit_df$OCDR_abs))^2)
R2 <- 1 - rss/tss
cat(sprintf("Pseudo R^2 on original scale: %.3f\n", R2))

# === 9) plot diagnostics ======================================================
# scatter (Jmag vs OCDR) and fitted curve (we'll use combined Jmag or predicted)
ggplot(fit_df, aes(x = Jmag_med, y = OCDR_abs)) +
  geom_point() +
  geom_line(aes(y = OCDR_pred), color = "red") +
  scale_y_continuous(name = "OCDR (m/s)") +
  scale_x_continuous(name = "J (dB)") +
  ggtitle("OCDR vs jiggle (power-dive seconds)")

# plot time series of OCDR and Jmag for inspection
ts <- merge(per_sec, fit_df[, .(second, OCDR_pred)], by="second", all.x=TRUE)
ggplot(ts, aes(x = second)) +
  geom_line(aes(y = OCDR_med), color="blue") +
  geom_line(aes(y = OCDR_pred), color="red", linetype=2) +
  geom_line(aes(y = (Jmag_med - median(Jmag_med, na.rm=TRUE))/10), color="darkgreen") + # scaled J for overlay
  ggtitle("Time series: OCDR (blue), predicted (red), scaled J (green)") +
  ylab("OCDR (m/s) / scaled J")

# Bland-Altman style check on log scale
plot(fit_df$logO, log(fit_df$OCDR_pred), xlab="log(OCDR_obs)", ylab="log(OCDR_pred)")
abline(0,1, col="red")

# Use absolute OCDR (remove negatives)
df <- power_secs[OCDR_abs > 0 & is.finite(Jmag_med)] %>%
  dplyr::filter(Jmag_med < -10)

# STARTING VALUES (important to avoid convergence issues)
a_start <- min(df$OCDR_abs, na.rm=TRUE)
b_start <- 0.05   # small positive slope is typical for dB vs speed

# === NONLINEAR EXPONENTIAL FIT (Cade et al. model) =============================
nls_fit <- nls(
  OCDR_abs ~ a * exp(b * Jmag_med),
  data = df,
  start = list(a = a_start, b = b_start)
)

summary(nls_fit)

# Extract coefficients
coef_a <- coef(nls_fit)[["a"]]
coef_b <- coef(nls_fit)[["b"]]

cat("Fitted model:\n")
cat(sprintf("  OCDR = %.4f * exp(%.4f * J)\n", coef_a, coef_b))

rss <- sum((df$OCDR_abs - df$OCDR_pred)^2)
tss <- sum((df$OCDR_abs - mean(df$OCDR_abs))^2)
R2 <- 1 - rss/tss
cat(sprintf("Pseudo R² = %.3f\n", R2))

# Predicted curve for plotting
df$OCDR_pred <- predict(nls_fit)

# === SCATTER PLOT IN CADE ET AL. STYLE =======================================
p <- ggplot(df, aes(x = Jmag_med, y = OCDR_abs)) +
  geom_point(alpha=0.6, size=2) +
  geom_line(aes(y = OCDR_pred), color="red", linewidth=1.2) +
  theme_classic(base_size = 14) +
  labs(
    x = "Jiggle (dB)",
    y = "OCDR speed (m/s)",
    title = "Jiggle vs Speed (Cade et al. style exponential fit)"
  ) +
  annotate(
    "text",
    x = min(df$Jmag_med) + 0.1 * diff(range(df$Jmag_med)),
    y = max(df$OCDR_abs) * 0.9,
    hjust = 0,
    label = sprintf("U = %.3f · e^[%.3f·J]", coef_a, coef_b),
    size = 5
    
  )

p

ggplot(df, aes(x = Jmag_med, y = OCDR_abs - OCDR_pred)) +
  geom_point() + geom_hline(yintercept=0, color="red") +
  theme_classic() +
  labs(x = "Jiggle (dB)", y = "Residual (m/s)")

# === 10) optional: robust fit / single-axis fits & comparisons ================
# single-axis fits:
lm_x <- lm(logO ~ Jx_med, data=fit_df); lm_y <- lm(logO ~ Jy_med, data=fit_df); lm_z <- lm(logO ~ Jz_med, data=fit_df)
summary(lm_x); summary(lm_y); summary(lm_z)

# compute predicted OCDR using only combined magnitude Jmag_med (one-var fit)
lm_mag <- lm(logO ~ Jmag_med, data = fit_df)
fit_df[, OCDR_pred_mag := exp(predict(lm_mag, newdata = fit_df))]

# compare R2 of multivariate vs magnitude-only
cat("R2 multivariate (original-scale):", R2, "\n")
rss_mag <- sum((fit_df$OCDR_abs - fit_df$OCDR_pred_mag)^2)
R2_mag <- 1 - rss_mag/tss
cat("R2 magnitude-only (original-scale):", R2_mag, "\n")

# === 11) Save results (optional) =============================================
per_sec$ID = o$ID[1]
fit_df$ID = o$ID[1]

# write per-second jiggle + OCDR summary for later analysis
fwrite(per_sec, file = paste0(output_path, "/", DT$ID[1], "_1Hz_jiggle.csv"))
fwrite(fit_df, file = paste0(output_path2, "/", DT$ID[1], "_jiggle_fit.csv"))

cat("Finished jiggle computation and preliminary regression.\n")
