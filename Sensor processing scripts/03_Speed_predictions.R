library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(broom)


setwd("~/SCIENCE/PhD/5. LittleLeo/Deployments")
# Set your top-level folder here

setwd("~/SCIENCE/PhD/5. LittleLeo/Deployments")
# Set your top-level folder here

data_path <- here::here("Processed", "02_Jiggle")

folder_path <- here::here("Processed", "02_Jiggle")

output_path <- here::here("Processed", "03_Speed")

# List all matching files recursively
files <- list.files(
  path = folder_path,
  pattern = "_jiggle\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

Deployment = 8

o <- read.csv(files[Deployment]) # cycle manually through all files


#a <- read.csv(paste0(data_path, "/LL240825/LL240825_DT_jiggle.csv"))
#b <- read.csv(paste0(data_path, "/LL250709/LL250709_DT_jiggle.csv"))
#c <- read.csv(paste0(data_path, "/LL250915/LL250915_DT_jiggle.csv"))

# ─────────────────────────────────────────────────────────────
# 0. PREP: Clean data structure
# ─────────────────────────────────────────────────────────────
# Expect columns:
# second, Jmag_med, Jx_med, Jy_med, Jz_med, OCDR_med, pitch_med, Depth_med, etc.

df <- o %>%
  rename(
    J = Jmag_med, 
    U = OCDR_med
  ) %>%
  mutate(
    logU = log(abs(U)),         # needed for linearised regression
    absU = abs(U),
    pitch_deg = pitch_med * 180/pi
  ) #%>%
  #dplyr::filter(is.finite(J), is.finite(U), absU > 0)


# ─────────────────────────────────────────────────────────────
# 1. OUTLIER REMOVAL
#    A: Remove non-power-dive bins |pitch| < 45° or NA
#    B: Remove extreme jiggle/speed outliers (IQR filtering)
# ─────────────────────────────────────────────────────────────

# A: keep only power-dive seconds
df_pow <- df %>%
  dplyr::filter(abs(pitch_deg) >= 45)

# B: remove IQR outliers in J and U
iqr_filter <- function(x, k=3) {  # 3×IQR filter (Cade-style robust)
  q <- quantile(x, probs=c(.25,.75), na.rm=TRUE)
  iqr <- q[2] - q[1]
  x > (q[1] - k*iqr) & x < (q[2] + k*iqr)
}

df_pow <- df_pow %>%
  dplyr::filter(iqr_filter(J), iqr_filter(absU))

# remove values from tag shedding for deployment LL240825
#df_pow <- df_pow %>%
#  dplyr::filter(second < 42869)

message("Remaining power-dive points: ", nrow(df_pow))


# ─────────────────────────────────────────────────────────────
# 2. FIT BOTH CADE REGRESSIONS
# ─────────────────────────────────────────────────────────────

# ------------------------------------------------------
# 2.1 Cade Model 1: Simple exponential
#     U = a * exp(b * J)
#     linearised: logU = log(a) + b*J
# ------------------------------------------------------

mod1 <- lm(logU ~ J, data=df_pow)
summary(mod1)

a1 <- exp(coef(mod1)[1])   # intercept in exponential form
b1 <- coef(mod1)[2]

# ------------------------------------------------------
# 2.2 Cade Model 2: Three-axis weighted version
#     U = a e^{ b(c1 Jx + c2 Jy + c3 Jz) }
#     with constraint: c1 + c2 + c3 = 1 (weights sum to 1)
#
#     We reparametrise:
#         W = c1*Jx + c2*Jy + (1-c1-c2)*Jz
#
#     and fit U = a e^{b W} as an nls model.
# ------------------------------------------------------

df_pow <- df_pow %>% 
  mutate(Jx = Jx_med, Jy = Jy_med, Jz = Jz_med)

# starting values
start_vals <- list(a = a1, b = b1, c1 = 0.33, c2 = 0.33)

mod2 <- nls(
  absU ~ a * exp( b * (c1*Jx + c2*Jy + (1 - c1 - c2)*Jz) ),
  data = df_pow,
  start = start_vals,
  control = nls.control(maxiter=200, warnOnly=TRUE)
)

summary(mod2)

coef_m2 <- coef(mod2)
a2 <- coef_m2["a"]
b2 <- coef_m2["b"]
c1 <- coef_m2["c1"]
c2 <- coef_m2["c2"]
c3 <- 1 - c1 - c2

# ─────────────────────────────────────────────────────────────
# 3. COMPARE REGRESSIONS
#    Using R² on linearised version for Model 1,
#    Using pseudo-R² for Model 2
# ─────────────────────────────────────────────────────────────

# R2 for model 1 (linear model)
R2_m1 <- summary(mod1)$r.squared

# Pseudo-R2 for model 2 (based on residuals)
df_pow <- df_pow %>%
  mutate(U_pred_m2 = predict(mod2))

rss2 <- sum((df_pow$absU - df_pow$U_pred_m2)^2)
tss2 <- sum((df_pow$absU - mean(df_pow$absU))^2)
R2_m2 <- 1 - rss2/tss2

comparison_table <- tibble(
  Model = c("Cade Model 1: U = a exp(bJ)",
            "Cade Model 2: U = a exp[b(c1Jx+c2Jy+c3Jz)]"),
  a = c(a1, a2),
  b = c(b1, b2),
  c1 = c(NA, c1),
  c2 = c(NA, c2),
  c3 = c(NA, c3),
  R2 = c(R2_m1, R2_m2)
)

comparison_table

# ─────────────────────────────────────────────────────────────
# 4. PREDICT SPEED FOR THE FULL DEPLOYMENT
#    Use the model with best R² (select automatically)
#    Provide predictions + 95% CI (model 1 only, linear model)
# ─────────────────────────────────────────────────────────────

#best_model <- if (R2_m1 >= R2_m2) "mod1" else "mod2"
#message("Best model is: ", best_model)
best_model <- "mod2"

df_pred <- df %>% mutate(U_pred = NA_real_, U_low = NA_real_, U_high = NA_real_)

if(best_model == "mod1") {
  
  # Predictions from linear model with confidence intervals
  preds <- predict(mod1, newdata=df, interval="confidence")
  
  df_pred <- df_pred %>%
    mutate(
      logU_pred = preds[,"fit"],
      logU_low  = preds[,"lwr"],
      logU_high = preds[,"upr"],
      U_pred    = exp(logU_pred),
      U_low     = exp(logU_low),
      U_high    = exp(logU_high)
    )
  
} else {
  
  # Nonlinear model — no built-in CI → use bootstrap or delta method
  # For now: point estimates only
  df_pred <- df_pred %>%
    mutate(
      W = c1*Jx_med + c2*Jy_med + c3*Jz_med,
      U_pred = a2 * exp(b2 * W)
    )
  
  message("Note: CI not computed for NLS model (need bootstrap).")
}

# Final predicted full time series
df_pred

# ─────────────────────────────────────────────────────────────
# 5. OPTIONAL PLOTTING
# ─────────────────────────────────────────────────────────────

# Create smooth prediction grid
J_seq <- seq(
  min(df_pow$J, na.rm = TRUE),
  max(df_pow$J, na.rm = TRUE),
  length.out = 200
)

pred_df <- data.frame(J = J_seq)

preds <- predict(
  mod1,
  newdata = pred_df,
  interval = "confidence",
  level = 0.95
)

pred_df <- pred_df %>%
  mutate(
    logU_fit  = preds[, "fit"],
    logU_lwr  = preds[, "lwr"],
    logU_upr  = preds[, "upr"],
    U_fit = exp(logU_fit),
    U_lwr = exp(logU_lwr),
    U_upr = exp(logU_upr)
  )

eq_label <- sprintf(
  "R² = %.2f",
  R2_m2
)

ggplot(df_pow, aes(J, absU)) +
  
  # raw data
  geom_point(
    alpha = 0.5,
    size = 2,
    colour = "black"
  ) +
  
  # confidence interval ribbon
  geom_ribbon(
    data = pred_df,
    aes(x = J, ymin = U_lwr, ymax = U_upr),
    fill = "grey70",
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  
  # fitted regression
  geom_line(
    data = pred_df,
    aes(x = J, y = U_fit),
    colour = "red",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  
  # equation + R² annotation
  annotate(
    "text",
    x = min(df_pow$J, na.rm = TRUE),
    y = max(df_pow$absU, na.rm = TRUE),
    label = eq_label,
    hjust = 0,
    vjust = 1,
    size = 6
  ) +
  
  theme_classic(base_size = 16) +
  
  labs(
    x = "Jiggle RMS amplitude (dB)",
    y = "OCDR (m/s)"
  )



ggplot(df_pow, aes(J, absU)) +
  geom_point(alpha=0.5) +
  geom_line(aes(y = predict(mod1, newdata=df_pow) %>% exp()), 
            col="red", size=1.1) +
  theme_classic() +
  labs(x="Jiggle RMS amplitude (dB)", y="OCDR (m/s)")

ggplot(df_pred, aes(second, U_pred)) +
  geom_ribbon(aes(x = second, ymin = U_low, ymax = U_high), color = "grey")+
  geom_line() +
  #xlim(0,32935)+
  theme_classic() +
  labs(title="Predicted speed for full deployment",
       x="Time (s)", y="Predicted speed (m/s)")

ggplot(df_pred, aes(U_pred))+
  geom_histogram(binwidth = 0.01, color = "black")+
  #xlim(c(0.4, 4))+
  theme_bw()+
  labs(title="Predicted speed for full deployment",
       x="Predicted speed (m/s)", y="Frequency")

summary(df_pred$U_pred)

write.csv(df_pred, paste0(output_path, "/", df$ID[1], "_speed.csv"), row.names = F)
write.csv(comparison_table, paste0(output_path, "/", df$ID[1], "_model_fit.csv"), row.names = F)

          