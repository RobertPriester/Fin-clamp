library(dplyr)
library(purrr)
library(readr)
library(lubridate)
library(hms)
library(drc)
library(runner)
library(ggplot2)
library(readxl)
#library(qqplotr)
#library(ggdist)

setwd("~/SCIENCE/PhD/5. LittleLeo/Deployments")
# Set your top-level folder here

# Defne paths
data_path <- here::here("Data")

merged_path <- here::here("Processed", "04_Merged")

result_path <- here::here("Results")

# Import data
all <- read.csv(paste0(merged_path, "/All_1Hz_merged_260223.csv"))

all <- all %>%
  mutate(
    Datetime = as.POSIXct(Datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    Year = as.integer(year(Datetime)),
    Deployment_start = as.POSIXct(Deployment_start, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    time_since_tagging_sec = as.numeric(Datetime - Deployment_start, units = "secs"),
    time_since_tagging_hms = as_hms(time_since_tagging_sec),
    Recovery_status = case_when(
      is.finite(time_since_tagging_sec) & time_since_tagging_sec > 7 * 60 * 60 ~ "Recovered",
      TRUE ~ "Not recovered")
  )

LL01 <- all %>%
  filter(ID == "LL240719")
LL02 <- all %>%
  filter(ID == "LL240825")
LL03 <- all %>%
  filter(ID == "LL240828")
LL04 <- all %>%
  filter(ID == "LL240920")
LL05 <- all %>%
  filter(ID == "LL250522")
LL06 <- all %>%
  filter(ID == "LL250709")
LL07 <- all %>%
  filter(ID == "LL250721")
LL08 <- all %>%
  filter(ID == "LL250915")

# Plot themes #####
pub_theme = theme_classic()+
  theme(#panel.grid.major = element_line(size = 0.5, color = "grey"),
    # axis.line = element_line(size = 0.7, color = "black"), 
    text = element_text(size = 14),
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 14),
    axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"))


# Functions ####

## Quantile help function
qq_sim_envelope <- function(x, nsim = 1000, conf = 0.99) {
  x <- na.omit(x)
  n <- length(x)
  
  z <- qnorm(ppoints(n))
  
  sims <- replicate(nsim, sort(rnorm(n, mean(x), sd(x))))
  
  alpha <- (1 - conf) / 2
  
  data.frame(
    theoretical = z,
    sample = sort(x),
    lower = apply(sims, 1, quantile, alpha),
    upper = apply(sims, 1, quantile, 1 - alpha)
  )
}

# General Overview #####

summary <- all %>%
  group_by(ID) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        min  = ~min(.x, na.rm = TRUE),
        mean = ~mean(.x, na.rm = TRUE),
        max  = ~max(.x, na.rm = TRUE)
      )
    )
  ) %>%
  ungroup()


# Recovery period #####
 
RecoveryDF <- all%>%
  mutate(time_passed = ((time_since_tagging_sec + 900) %/% 900) * 900) %>%
  mutate(time_passed = time_passed/60/60) %>%
  #filter(Surface_X == 0) %>%
  #filter(Depth > 2) %>%
  filter(DominantCycle > 0.6) %>%
  group_by(ID, time_passed) %>%
  #summarise(TBC_X = mean(DominantCycle_X, na.rm = T),
  #          TBC_Y = mean(DominantCycle_Y, na.rm = T),
  #          TBC_Z = mean(DominantCycle_Z, na.rm = T),
  #          count = n()) %>%
  summarise(TBC_Z = mean(DominantCycle, na.rm = T),
            count = n()) %>%
  mutate(time_passed = as.numeric(time_passed)) %>% 
  ungroup() %>%
  filter(count > 100)

ggplot(RecoveryDF)+
  geom_point(aes(time_passed, TBC_Z, color = ID))

## 01 ####

Recovery_DF <- LL01 %>%
  mutate(time_passed = ((time_since_tagging_sec + 900) %/% 900) * 900) %>%
  mutate(time_passed = time_passed/60/60) %>%
  #filter(Surface_X == 0) %>%
  #filter(Depth_3s > 2) %>%
  group_by(time_passed) %>%
  #summarise(TBC_X = mean(DominantCycle_X, na.rm = T),
  #          TBC_Y = mean(DominantCycle_Y, na.rm = T),
  #          TBC_Z = mean(DominantCycle_Z, na.rm = T),
  #          count = n()) %>%
  summarise(TBC_Z = mean(DominantCycle, na.rm = T),
            count = n()) %>%
  mutate(time_passed = as.numeric(time_passed)) %>% 
  ungroup() %>%
  filter(count > 100)


model <- drm(TBC_Z ~ time_passed, fct = AR.3(names = c("init", "plateau", "m")), data = Recovery_DF)

summary(model)

plot(model, log = '', main = "Asymptotic regression")

# Extract fitted values and residuals
fitted_values <- fitted(model)
residuals <- residuals(model, type = 'standardised')

# Plot residuals
ggplot(data = NULL, aes(x = fitted_values, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  ylab('Residuals') + xlab('Fitted Values') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Create a dataframe with residuals
residuals_df <- data.frame(residuals = residuals)


qq_df <- qq_sim_envelope(residuals_df$residuals)

ggplot(qq_df, aes(theoretical, sample)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "grey70", alpha = 0.5) +
  geom_abline(intercept = mean(residuals_df$residuals),
              slope = sd(residuals_df$residuals)) +
  geom_point(size = 1) +
  labs(x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_bw()


# Check autocorrelation
acf(residuals, main = "Residual Autocorrelation")

# Extract Plateau (right asymptote)
Plateau <- as.vector(model$coefficients['plateau:(Intercept)'])

# Calculate the tail beat cycle where it reached 80% of the difference between the initial 
# post-release value and the fully recovered value

RecoverdCycle <- 
  Recovery_DF$TBC_Z[which(Recovery_DF$time_passed == 0.25)] + 
  (Plateau - Recovery_DF$TBC_Z[which(Recovery_DF$time_passed == 0.25)])*0.8

# Find the point in time
RecoveryTime <- uniroot(function(x) predict(model, data.frame(time_passed = x)) - RecoverdCycle,
                        interval = c(0, max(Recovery_DF$time_passed)))$root


# Plotting

pm <- predict(model, data.frame(time_passed = seq(0, (max(Recovery_DF$time_passed) + 0.25), 0.25)),
              interval="confidence")

pldt <- data.frame(time_passed = seq(0, (max(Recovery_DF$time_passed) + 0.25), 0.25)) %>%
  mutate(p = pm[,1], pmin = pm[,2], pmax = pm[,3])


p2 <- ggplot() +
  geom_point(data = Recovery_DF, aes(x = time_passed, y = TBC_Z), size = 2) +
  geom_ribbon(data = pldt, aes(x=time_passed, y = p, ymin = pmin, ymax = pmax), alpha=0.2) +
  geom_line(data = pldt, aes(x = time_passed, y = p), size = 1, colour = "red") + 
  ylab('Tailbeat Cycle (sec)') + xlab('Time after tagging (hours)') +
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) +
  ggtitle('#1 (F, 125 cm TL, Rod-and-Reel capture)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title =  element_text(size = 14),
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(vjust = 2, size = 12), 
        axis.ticks.length = unit(2,'mm'),
        axis.text = element_text(size = 10, colour = 'black'), 
        axis.ticks = element_line(colour = "black", linewidth = 0.25),
        plot.margin = margin(2.5, 2.5, 2.5 , 2.5, "mm"),
        axis.line = element_line(colour = "black", 
                                 linewidth = 0.25, linetype = "solid")) +
  geom_vline(xintercept = RecoveryTime, size = 1, colour = 'red', linetype = 'dashed') +
  geom_text(aes(x=RecoveryTime, y = 0.5, label= paste0('80% Asymptote +/- ', round(RecoveryTime, digits = 2), ' hours')), 
            size = 4, angle=90, vjust=-0.5, hjust=0) +
  pub_theme

p2

lin_mod <- lm(TBC_Z ~ time_passed, data = Recovery_DF)

AIC(model)
AIC(lin_mod)

## 02 ####
Recovery_DF <- LL02 %>%
  mutate(time_passed = ((time_since_tagging_sec + 900) %/% 900) * 900) %>%
  mutate(time_passed = time_passed/60/60) %>%
  #filter(Surface_X == 0) %>%
  #filter(Depth_3s > 2) %>%
  group_by(time_passed) %>%
  #summarise(TBC_X = mean(DominantCycle_X, na.rm = T),
  #          TBC_Y = mean(DominantCycle_Y, na.rm = T),
  #          TBC_Z = mean(DominantCycle_Z, na.rm = T),
  #          count = n()) %>%
  summarise(TBC_Z = mean(DominantCycle, na.rm = T),
            count = n()) %>%
  mutate(time_passed = as.numeric(time_passed)) %>% 
  ungroup() %>%
  filter(count > 100)


model <- drm(TBC_Z ~ time_passed, fct = AR.3(names = c("init", "plateau", "m")), data = Recovery_DF)

summary(model)

lin_mod <- lm(TBC_Z ~ time_passed, data = Recovery_DF)

AIC(model)
AIC(lin_mod)

plot(model, log = '', main = "Asymptotic regression")

# Extract fitted values and residuals
fitted_values <- fitted(model)
residuals <- residuals(model, type = 'standardised')

# Plot residuals
ggplot(data = NULL, aes(x = fitted_values, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  ylab('Residuals') + xlab('Fitted Values') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Create a dataframe with residuals
residuals_df <- data.frame(residuals = residuals)


qq_df <- qq_sim_envelope(residuals_df$residuals)

ggplot(qq_df, aes(theoretical, sample)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "grey70", alpha = 0.5) +
  geom_abline(intercept = mean(residuals_df$residuals),
              slope = sd(residuals_df$residuals)) +
  geom_point(size = 1) +
  labs(x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_bw()


# Check autocorrelation
acf(residuals, main = "Residual Autocorrelation")

# Extract Plateau (right asymptote)
Plateau <- as.vector(model$coefficients['plateau:(Intercept)'])

# Calculate the tail beat cycle where it reached 80% of the difference between the initial 
# post-release value and the fully recovered value

RecoverdCycle <- 
  Recovery_DF$TBC_Z[which(Recovery_DF$time_passed == 0.25)] + 
  (Plateau - Recovery_DF$TBC_Z[which(Recovery_DF$time_passed == 0.25)])*0.8

# Find the point in time
RecoveryTime <- uniroot(function(x) predict(model, data.frame(time_passed = x)) - RecoverdCycle,
                        interval = c(0, max(Recovery_DF$time_passed)))$root

# In case above function returns error - run this code to see if asymptote was ever reached
f <- function(x)
  predict(model, data.frame(time_passed = x)) - RecoverdCycle

f(0)
f(max(Recovery_DF$time_passed))

tp <- seq(0, max(Recovery_DF$time_passed), length.out = 1000)
pred <- predict(model, data.frame(time_passed = tp))

RecoveryTime <- tp[which.min(abs(pred - RecoverdCycle))]

# Plotting

pm <- predict(model, data.frame(time_passed = seq(0, (max(Recovery_DF$time_passed) + 0.25), 0.25)),
              interval="confidence")

pldt <- data.frame(time_passed = seq(0, (max(Recovery_DF$time_passed) + 0.25), 0.25)) %>%
  mutate(p = pm[,1], pmin = pm[,2], pmax = pm[,3])


p2 <- ggplot() +
  geom_point(data = Recovery_DF, aes(x = time_passed, y = TBC_Z), size = 2) +
  geom_ribbon(data = pldt, aes(x=time_passed, y = p, ymin = pmin, ymax = pmax), alpha=0.2) +
  geom_line(data = pldt, aes(x = time_passed, y = p), size = 1, colour = "red") + 
  ylab('Tailbeat Cycle (sec)') + xlab('Time after tagging (hours)') +
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) +
  ggtitle('#2 (F, 114 cm TL, Longline capture)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title =  element_text(size = 14),
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(vjust = 2, size = 12), 
        axis.ticks.length = unit(2,'mm'),
        axis.text = element_text(size = 10, colour = 'black'), 
        axis.ticks = element_line(colour = "black", linewidth = 0.25),
        plot.margin = margin(2.5, 2.5, 2.5 , 2.5, "mm"),
        axis.line = element_line(colour = "black", 
                                 linewidth = 0.25, linetype = "solid")) +
  #geom_vline(xintercept = RecoveryTime, size = 1, colour = 'red', linetype = 'dashed') +
  #geom_text(aes(x=RecoveryTime, y = 0.5, label= paste0('80% Asymptote +/- ', round(RecoveryTime, digits = 2), ' hours')), 
  #          size = 4, angle=90, vjust=-0.5, hjust=0) +
  pub_theme

p2

## 03 #####

Recovery_DF <- LL03 %>%
  mutate(time_passed = ((time_since_tagging_sec + 900) %/% 900) * 900) %>%
  mutate(time_passed = time_passed/60/60) %>%
  #filter(Surface_X == 0) %>%
  #filter(Depth_3s > 2) %>%
  group_by(time_passed) %>%
  #summarise(TBC_X = mean(DominantCycle_X, na.rm = T),
  #          TBC_Y = mean(DominantCycle_Y, na.rm = T),
  #          TBC_Z = mean(DominantCycle_Z, na.rm = T),
  #          count = n()) %>%
  summarise(TBC_Z = mean(DominantCycle, na.rm = T),
            count = n()) %>%
  mutate(time_passed = as.numeric(time_passed)) %>% 
  ungroup() %>%
  filter(count > 100)


model <- drm(TBC_Z ~ time_passed, fct = AR.3(names = c("init", "plateau", "m")), data = Recovery_DF)

summary(model)

lin_mod <- lm(TBC_Z ~ time_passed, data = Recovery_DF)

AIC(model)
AIC(lin_mod)

plot(model, log = '', main = "Asymptotic regression")

# Extract fitted values and residuals
fitted_values <- fitted(model)
residuals <- residuals(model, type = 'standardised')

# Plot residuals
ggplot(data = NULL, aes(x = fitted_values, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  ylab('Residuals') + xlab('Fitted Values') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Create a dataframe with residuals
residuals_df <- data.frame(residuals = residuals)


qq_df <- qq_sim_envelope(residuals_df$residuals)

ggplot(qq_df, aes(theoretical, sample)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "grey70", alpha = 0.5) +
  geom_abline(intercept = mean(residuals_df$residuals),
              slope = sd(residuals_df$residuals)) +
  geom_point(size = 1) +
  labs(x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_bw()


# Check autocorrelation
acf(residuals, main = "Residual Autocorrelation")

# Extract Plateau (right asymptote)
Plateau <- as.vector(model$coefficients['plateau:(Intercept)'])

# Calculate the tail beat cycle where it reached 80% of the difference between the initial 
# post-release value and the fully recovered value

RecoverdCycle <- 
  Recovery_DF$TBC_Z[which(Recovery_DF$time_passed == 0.25)] + 
  (Plateau - Recovery_DF$TBC_Z[which(Recovery_DF$time_passed == 0.25)])*0.8

# Find the point in time
RecoveryTime <- uniroot(function(x) predict(model, data.frame(time_passed = x)) - RecoverdCycle,
                        interval = c(0, max(Recovery_DF$time_passed)))$root


# Plotting

pm <- predict(model, data.frame(time_passed = seq(0, (max(Recovery_DF$time_passed) + 0.25), 0.25)),
              interval="confidence")

pldt <- data.frame(time_passed = seq(0, (max(Recovery_DF$time_passed) + 0.25), 0.25)) %>%
  mutate(p = pm[,1], pmin = pm[,2], pmax = pm[,3])


p2 <- ggplot() +
  geom_point(data = Recovery_DF, aes(x = time_passed, y = TBC_Z), size = 2) +
  geom_ribbon(data = pldt, aes(x=time_passed, y = p, ymin = pmin, ymax = pmax), alpha=0.2) +
  geom_line(data = pldt, aes(x = time_passed, y = p), size = 1, colour = "red") + 
  ylab('Tailbeat Cycle (sec)') + xlab('Time after tagging (hours)') +
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) +
  ggtitle('#3 (M, 97 cm TL, Longline capture)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title =  element_text(size = 14),
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(vjust = 2, size = 12), 
        axis.ticks.length = unit(2,'mm'),
        axis.text = element_text(size = 16, colour = 'black'), 
        axis.ticks = element_line(colour = "black", linewidth = 0.25),
        plot.margin = margin(2.5, 2.5, 2.5 , 2.5, "mm"),
        axis.line = element_line(colour = "black", 
                                 linewidth = 0.25, linetype = "solid")) +
  geom_vline(xintercept = RecoveryTime, size = 1, colour = 'red', linetype = 'dashed') +
  geom_text(aes(x=RecoveryTime, y = 0.5, label= paste0('80% Asymptote +/- ', round(RecoveryTime, digits = 2), ' hours')), 
            size = 5, angle=90, vjust=-0.5, hjust=0) +
  pub_theme

p2



## 07 ####
Recovery_DF <- LL07 %>%
  mutate(time_passed = ((time_since_tagging_sec + 900) %/% 900) * 900) %>%
  mutate(time_passed = time_passed/60/60) %>%
  #filter(Surface_X == 0) %>%
  #filter(Depth_3s > 2) %>%
  group_by(time_passed) %>%
  #summarise(TBC_X = mean(DominantCycle_X, na.rm = T),
  #          TBC_Y = mean(DominantCycle_Y, na.rm = T),
  #          TBC_Z = mean(DominantCycle_Z, na.rm = T),
  #          count = n()) %>%
  summarise(TBC_Z = mean(DominantCycle, na.rm = T),
            count = n()) %>%
  mutate(time_passed = as.numeric(time_passed)) %>% 
  ungroup() %>%
  filter(count > 100)


model <- drm(TBC_Z ~ time_passed, fct = AR.3(names = c("init", "plateau", "m")), data = Recovery_DF)

summary(model)

lin_mod <- lm(TBC_Z ~ time_passed, data = Recovery_DF)

AIC(model)
AIC(lin_mod)

plot(model, log = '', main = "Asymptotic regression")

# Extract fitted values and residuals
fitted_values <- fitted(model)
residuals <- residuals(model, type = 'standardised')

# Plot residuals
ggplot(data = NULL, aes(x = fitted_values, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  ylab('Residuals') + xlab('Fitted Values') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Create a dataframe with residuals
residuals_df <- data.frame(residuals = residuals)


qq_df <- qq_sim_envelope(residuals_df$residuals)

ggplot(qq_df, aes(theoretical, sample)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "grey70", alpha = 0.5) +
  geom_abline(intercept = mean(residuals_df$residuals),
              slope = sd(residuals_df$residuals)) +
  geom_point(size = 1) +
  labs(x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_bw()


# Check autocorrelation
acf(residuals, main = "Residual Autocorrelation")

# Extract Plateau (right asymptote)
Plateau <- as.vector(model$coefficients['plateau:(Intercept)'])

# Calculate the tail beat cycle where it reached 80% of the difference between the initial 
# post-release value and the fully recovered value

RecoverdCycle <- 
  Recovery_DF$TBC_Z[which(Recovery_DF$time_passed == 0.25)] + 
  (Plateau - Recovery_DF$TBC_Z[which(Recovery_DF$time_passed == 0.25)])*0.8

# Find the point in time
RecoveryTime <- uniroot(function(x) predict(model, data.frame(time_passed = x)) - RecoverdCycle,
                        interval = c(0, max(Recovery_DF$time_passed)))$root


# Plotting

pm <- predict(model, data.frame(time_passed = seq(0, (max(Recovery_DF$time_passed) + 0.25), 0.25)),
              interval="confidence")

pldt <- data.frame(time_passed = seq(0, (max(Recovery_DF$time_passed) + 0.25), 0.25)) %>%
  mutate(p = pm[,1], pmin = pm[,2], pmax = pm[,3])


p2 <- ggplot() +
  geom_point(data = Recovery_DF, aes(x = time_passed, y = TBC_Z), size = 2) +
  geom_ribbon(data = pldt, aes(x=time_passed, y = p, ymin = pmin, ymax = pmax), alpha=0.2) +
  geom_line(data = pldt, aes(x = time_passed, y = p), size = 1, colour = "red") + 
  ylab('Tailbeat Cycle (sec)') + xlab('Time after tagging (hours)') +
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) +
  ggtitle('#7 (M, 125.5 cm TL, Longline capture)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title =  element_text(size = 14),
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(vjust = 2, size = 12), 
        axis.ticks.length = unit(2,'mm'),
        axis.text = element_text(size = 10, colour = 'black'), 
        axis.ticks = element_line(colour = "black", linewidth = 0.25),
        plot.margin = margin(2.5, 2.5, 2.5 , 2.5, "mm"),
        axis.line = element_line(colour = "black", 
                                 linewidth = 0.25, linetype = "solid")) +
  geom_vline(xintercept = RecoveryTime, size = 1, colour = 'red', linetype = 'dashed') +
  geom_text(aes(x=RecoveryTime, y = 0.5, label= paste0('80% Asymptote +/- ', round(RecoveryTime, digits = 2), ' hours')), 
            size = 4.5, angle=90, vjust=-0.5, hjust=0) +
  pub_theme

p2


## ALL #####

# Function to compute metrics for one deployment
compute_metrics <- function(data, recoverd_cycle) {
  mod <- drm(TBC_Z ~ time_passed,
             fct = AR.3(names = c("init", "plateau", "m")),
             data = data)
  
  y <- data$TBC_Z
  y_hat <- fitted(mod)
  
  RSS <- sum((y - y_hat)^2)
  TSS <- sum((y - mean(y))^2)
  R2 <- 1 - RSS / TSS
  RMSE_val <- sqrt(mean((y - y_hat)^2))
  
  t_end <- max(data$time_passed)
  pred_end <- predict(mod, data.frame(time_passed = t_end))
  
  # Correctly extract plateau (asymptote)
  plateau <- coef(mod)[2]
  asymptote_dev <- abs(plateau - pred_end)
  
  # Robust recovery time
  tp <- seq(0, t_end, length.out = 1000)
  pred_tp <- predict(mod, data.frame(time_passed = tp))
  
  if (all(pred_tp < recoverd_cycle)) {
    RecoveryTime <- NA
  } else {
    idx <- which(diff(sign(pred_tp - recoverd_cycle)) != 0)
    if (length(idx) == 0) {
      RecoveryTime <- NA
    } else {
      RecoveryTime <- uniroot(
        function(x) predict(mod, data.frame(time_passed = x)) - recoverd_cycle,
        interval = c(tp[idx[1]], tp[idx[1] + 1])
      )$root
    }
  }
  
  Asymptote <- plateau
  
  tibble(R2 = R2,
         RMSE = RMSE_val,
         asymptote_dev = asymptote_dev,
         RecoveryTime = RecoveryTime,
         Asymptote = Asymptote)
}



Recovery_DF_all <- all %>%
  mutate(time_passed = ((time_since_tagging_sec + 900) %/% 900) * 900) %>%
  mutate(time_passed = time_passed/60/60) %>%
  #filter(Surface_X == 0) %>%
  #filter(Depth_3s > 2) %>%
  group_by(ID, time_passed) %>%
  #summarise(TBC_X = mean(DominantCycle_X, na.rm = T),
  #          TBC_Y = mean(DominantCycle_Y, na.rm = T),
  #          TBC_Z = mean(DominantCycle_Z, na.rm = T),
  #          count = n()) %>%
  summarise(TBC_Z = mean(DominantCycle, na.rm = T),
            count = n()) %>%
  mutate(time_passed = as.numeric(time_passed)) %>% 
  ungroup() %>%
  filter(count > 100)

# Apply to all deployments (assuming column 'DeploymentID')
summary_metrics <- Recovery_DF_all %>%
  group_by(ID) %>%
  group_modify(~compute_metrics(.x, RecoverdCycle))





## Dive variance ####

LL01 <- combined_tbf_d %>%
  filter(DeploymentID == "LL240825")

LL01_rec <- LL01 %>%
  mutate(time_since_tagging_hour = time_since_tagging_sec/60/60) %>%
  group_by(floor(time_since_tagging_hour)) %>%
  summarize(
    dive_var = var(Depth)
  )

# Tag shift ####

ggplot(all)+
  geom_line(aes(time_since_tagging_hms, pitch_offset_gam*180/pi, group = ID, color = Year))

ggplot(all)+
  geom_point(aes(time_since_tagging_sec/60/60, pitch_vv0*180/pi, fill = ID),shape = 21, alpha = 0.005)+
  geom_line(aes(time_since_tagging_sec/60/60, pitch_offset_gam*180/pi, group = ID, color = ID), linewidth = 0.8, alpha = 0.9)+
  #scale_color_viridis_d(option="turbo")+
  scale_color_manual(values = c("orange", "orange1", "orange2", "orange3",
                                "royalblue", "royalblue1", "royalblue2", "royalblue3"))+
  ylim(c(-18, 5))+
  labs(
    x = "Time since release (h)",
    y = expression("Pitch offset (" * degree * ")")
  ) +
  pub_theme+
  theme(
    legend.position = "none")

ggplot(tagshift_15min)+
  geom_point(aes(time_since_tagging_sec/60/60, abs(roll*180/pi), group = ID, color = ID), alpha = 0.2)+
  geom_smooth(aes(time_since_tagging_sec/60/60, abs(roll*180/pi), group = ID, color = ID), method = "gam", alpha = 0.9, se = F)+
  scale_color_manual(values = c("orange", "orange1", "orange2", "orange3",
                                "royalblue", "royalblue1", "royalblue2", "royalblue3"))+
  #ylim(c(-18, 5))+
  labs(
    x = "Time since release (h)",
    y = expression("Roll (" * degree * ")")
  ) +
  pub_theme+
  theme(
    legend.position = "none")


ggplot(LL02)+
  geom_line(aes(time_since_tagging_hms, pitch_offset_gam*180/pi, group = ID), linewidth = 1)+
  geom_point(aes(time_since_tagging_hms, pitch_vv0*180/pi), alpha = 0.05)+
  geom_line(aes(time_since_tagging_hms, Jerk, group = ID))

ggplot(LL03)+
  geom_point(aes(time_since_tagging_hms, pitch_vv0*180/pi), alpha = 0.05)+
  geom_line(aes(time_since_tagging_hms, VeDBA, color = VeDBA))+
  geom_line(aes(time_since_tagging_hms, pitch_offset_gam*180/pi), linewidth = 1.2)+
  scale_color_viridis_c(option="turbo")+
  labs(
    x = "Time since release (h)",
    y = "Pitch offset (degrees)"
  ) +
  pub_theme+
  theme(
    legend.position = "none")
  #geom_line(aes(time_since_tagging_hms, DominantCycleAmplitude*100))

tagshift_15min <- all %>%
  mutate(time_since_tagging_15min = floor(time_since_tagging_sec/60/15)) %>%
  group_by(ID, time_since_tagging_15min, TL) %>%
  summarize(
    roll = mean(roll),
    pitch = mean(pitch_offset_gam),
    jerk = max(Jerk, na.rm = F),
    time_since_tagging_sec = mean(time_since_tagging_sec)
  ) %>%
  ungroup()



ggplot(tagshift_15min)+
  geom_point(aes(time_since_tagging_hms, abs(pitch*180/pi)+50, group = ID, color = jerk))+
  geom_smooth(aes(time_since_tagging_hms, abs(pitch*180/pi)+50, group = ID, color = jerk), method = "gam")

tagshift <- all %>%tagshift <- all %>%tagshift <- all %>%
  group_by(ID) %>%
  arrange(ID, RunNo) %>%   # ensure correct order
  mutate(pitch_offset_rate = abs(pitch_offset_gam - lag(pitch_offset_gam))) %>%
  ungroup()

ggplot(tagshift)+
  geom_line(aes(time_since_tagging_hms, pitch_offset_rate*180/pi*60*24, group = ID, color = ID))

ggplot(tagshift)+
  geom_point(aes(Jerk, pitch_offset_rate*180/pi), alpha = 0.05)

## Summary Table ####

tagshift_summary <- all %>%
  filter(!is.na(pitch_offset_gam)) %>%
  group_by(ID) %>%
  summarize(
    pitch_start = first(pitch_offset_gam*180/pi, na.rm = T),
    pitch_end = last(pitch_offset_gam*180/pi, na.rm = T),
    pitch_min = min(pitch_offset_gam*180/pi, na.rm = T),
    pitch_max = max(pitch_offset_gam*180/pi, na.rm = T),
    pitch_mean = mean(pitch_offset_gam*180/pi, na.rm = T),
    pitch_abs_diff = pitch_max - pitch_min,
    pitch_start_end = pitch_end - pitch_start,
    pitch_shift_per_hour = pitch_start_end/(max(second)/60/60))

write.csv(tagshift_summary, paste0(result_path, "/Tagshift/Tagshift_summary.csv"), row.names = F)

# TBF ~ HBF #####

tbf <- all %>%
  dplyr::select(ID, Year, Depth_3s, DominantHz, TBF_annotated) %>%
  filter(!is.na(Year) & TBF_annotated < 3) %>%
  filter(between(TBF_annotated, 0.6, 1.6)) %>% # remove all values outside of CWT range
  mutate(class = case_when(
    Depth_3s < 2 ~ "surface",
    TRUE      ~ "deep"
  ))

sum <- tbf %>%
  group_by(ID, class) %>%
  summarize(
    n = n()
  )


### Prototype vs. Miniaturized ####
lm_stats <- tbf %>%
  group_by(Year) %>%
  do({
    model <- lm(TBF_annotated ~ DominantHz, data = .)  # <-- axes switched here
    data.frame(
      r  = cor(.$DominantHz, .$TBF_annotated),
      r2 = summary(model)$r.squared,
      slope = coef(model)[2],
      intercept = coef(model)[1]
    )
  })


ggplot(tbf, aes(x = TBF_annotated, y = DominantHz))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8)+
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~ Year, scales = "fixed",
             labeller = labeller(
               Year = c(
                 "2024" = "Prototype",
                 "2025" = "Miniaturized"
               )
             )) +
  geom_text(
    data = lm_stats,
    aes(
      label = paste0(
        "r = ", round(r, 2), "\n",
        "R² = ", round(r2, 2), "\n",
        "y = ", round(slope, 2), "x + ", round(intercept, 2)
      )
    ),
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.3,
    inherit.aes = FALSE,
    size = 4
  ) +
  labs(
    x = "Visual HBF",
    y = "Accelerometry TBF")+
  pub_theme


### Surface vs. Deep #####
lm_stats <- tbf %>%
  group_by(class) %>%
  do({
    model <- lm(TBF_annotated ~ DominantHz, data = .)  # <-- axes switched here
    data.frame(
      r  = cor(.$DominantHz, .$TBF_annotated),
      r2 = summary(model)$r.squared,
      slope = coef(model)[2],
      intercept = coef(model)[1]
    )
  })


ggplot(tbf, aes(x = TBF_annotated, y = DominantHz))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8)+
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~ class, scales = "fixed",
             labeller = labeller(
               class = c(
                 "deep" = "a - Deeper than 2 m",
                 "surface" = "b - Shallower than 2 m"
               )
             )) +
  geom_text(
    data = lm_stats,
    aes(
      label = paste0(
        "r = ", round(r, 2), "\n",
        "R² = ", round(r2, 2), "\n",
        "y = ", round(slope, 2), "x + ", round(intercept, 2)
      )
    ),
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.3,
    inherit.aes = FALSE,
    size = 4
  ) +
  labs(
    x = "Visual HBF",
    y = "Accelerometry TBF")+
  pub_theme+
  theme(strip.text = element_text(size = 16))

### All deployments #####
lm_stats <- tbf %>%
  dplyr::filter(ID != "LL240719") %>%
  group_by(ID) %>%
  do({
    model <- lm(TBF_annotated ~ DominantHz, data = .)  # <-- axes switched here
    data.frame(
      r  = cor(.$DominantHz, .$TBF_annotated),
      r2 = summary(model)$r.squared,
      slope = coef(model)[2],
      intercept = coef(model)[1]
    )
  })


ggplot(tbf %>%
         dplyr::filter(ID != "LL240719")
       , aes(x = TBF_annotated, y = DominantHz))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8)+
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", se = FALSE)+
  facet_wrap(~ ID, scales = "fixed",
             labeller = labeller(
               ID = c(
                 "LL240825" = "#2",
                 "LL240920" = "#4",
                 "LL250522" = "#5",
                 "LL250721" = "#7"
               )
             )) +
  geom_text(
    data = lm_stats,
    aes(
      label = paste0(
        "r = ", round(r, 2), "\n",
        "R² = ", round(r2, 2), "\n",
        "y = ", round(slope, 2), "x + ", round(intercept, 2)
      )
    ),
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.3,
    inherit.aes = FALSE,
    size = 4
  ) +
  labs(
    x = "Visual HBF",
    y = "Accelerometry TBF")+
  pub_theme+
  theme(strip.text = element_text(size = 16))


## Depth profile ####

ggplot(LL05)+
  geom_line(aes(Datetime, Depth*-1, color = Temp), linewidth = 1.1)+
  scale_color_viridis_c(option="plasma", name = "Temp. (°C)")+
  xlim(c(min(LL05$Datetime), max(LL05$Datetime)))+
  scale_x_datetime(expand = c(0, 0)) +
  labs(
    x = "Time",
    y = "Depth (m)")+
  pub_theme#+
  theme(
    legend.position = "none")


## TBF quality #####

LL07 <- combined_tbf_d %>%
  filter(DeploymentID == "LL250721")

ggplot(LL07, aes(Depth, EntropyOfFrequency))+
  geom_boxplot(aes(group = Depth), alpha = 0.05)


## Speed analysis ####

ggplot(all)+
  geom_histogram(aes(U_pred, fill = ID), binwidth = 0.01)+
  xlim(0.1, 3)

speed_summary <- all %>%
  group_by(ID) %>%
  summarize(
    n = n(),
    min = min(U_pred, na.rm = T),
    mean = mean(U_pred, na.rm = T),
    median = median(U_pred, na.rm = T),
    max = max(U_pred, na.rm = T),
    perc02.5 = quantile(U_pred, probs = 0.025, na.rm = T),
    perc95 = quantile(U_pred, probs = 0.95, na.rm = T),
    perc97.5 = quantile(U_pred, probs = 0.975, na.rm = T),
    perc99 = quantile(U_pred, probs = 0.99, na.rm = T),
    prob_above2 = mean(U_pred >1.5, na.rm = T)
    )



combined_jiggle <- combined_jiggle %>%
  left_join(speed_summary %>% dplyr::select(DeploymentID, perc95), by = "DeploymentID")

burst <- combined_jiggle %>%
  filter(U_pred > perc95)

burst_summary <- burst %>%
  group_by(DeploymentID) %>%
  summarise(
    n_burst = n()) %>%
  right_join(speed_summary) %>%
  mutate(n_burst/n*100)



ggplot(burst)+
  geom_histogram(aes(U_pred, fill = DeploymentID), binwidth = 0.01)+
  xlim(0.1, 3)
