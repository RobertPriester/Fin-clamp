library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)



# Load  files
sims <- read.csv("AllSimulations.csv")

meta <- read_excel("Result_summary.xlsx")

mesh <- read_excel("Result_summary.xlsx", sheet = "Mesh dependancy")

### Plot themes #####
pub_theme = theme_classic()+
  theme(#panel.grid.major = element_line(size = 0.5, color = "grey"),
    # axis.line = element_line(size = 0.7, color = "black"), 
    text = element_text(size = 14),
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 14),
    axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"))

cfd <- sims %>%
  filter(Time == 2000) %>%
  mutate(CD = TOTAL_FORCE_X/(0.5*1024.065*(Velocity^2)*Ref_area),
         Power = TOTAL_FORCE_X*Velocity)



cfd_clean <- cfd %>%
  select(Shark_size, Velocity, Tag, TOTAL_FORCE_X, CD, Power)

cfd_wide <- cfd_clean %>%
  pivot_wider(
    names_from = Tag,
    values_from = c(TOTAL_FORCE_X, CD, Power)
  ) %>%
  mutate(
    delta_drag_N = TOTAL_FORCE_X_YES - TOTAL_FORCE_X_NO,
    delta_Cd     = CD_YES - CD_NO,
    rel_drag_pct = delta_drag_N / TOTAL_FORCE_X_NO * 100,
    delta_power = Power_YES - Power_NO,
    rel_power_pct = delta_power/Power_NO*100
  ) %>%
  left_join(cfd %>%
              filter(Tag == "YES") %>%
              select(Shark_size, Velocity, Tag.Shark_weight, Tag.shark_area), by = c("Shark_size", "Velocity")) %>%
  mutate(
    Tag.Shark_weight = Tag.Shark_weight*100,
    Tag.shark_area = Tag.shark_area*100
  )

#write.csv(cfd_wide, "C:/Users/Lenovo/Documents/SCIENCE/PhD/5. LittleLeo/Drag modelling/04_results/CFD_summary.csv", row.names = F)

cfd_wide_1 <- cfd_wide %>%
  filter(Velocity == 1)

sims <- sims %>%
  group_by(Run) %>%
  arrange(Time, .by_group = TRUE) %>%
  mutate(
    CD_diff = FORCE_COEFFICIENT_CD - last(FORCE_COEFFICIENT_CD)
  ) %>%
  ungroup()


mesh <- mesh %>%
  filter(
    Run == "Manual"
  ) %>%
  mutate(
    CD = TOTAL_FORCE_X/(0.5*1024.065*0.00963700),
    CD_diff = CD-first(CD)
  )

# Simulation stability #####
ggplot(sims, aes(Time, abs(CD_diff)))+
  geom_line(aes(group = Run), alpha = 0.15)+
  coord_cartesian(xlim = c(0, 2000),
                  ylim = c(0, 0.2))+
  geom_hline(yintercept = 0.01, color = "red", linetype = "dashed", alpha = 0.6, linewidth = 1)+
  labs(x = "Time (s)",
       y = expression("Absolute change in drag coefficient (" * abs(Delta * C[D]) * ")")) +
  pub_theme
  
# Mesh independancy ####
ggplot(mesh, aes(Volumes, CD_diff))+
  geom_line(linewidth = 1, alpha = 0.8)+
  geom_point(size = 2, alpha = 0.8)+
  geom_hline(yintercept = 0.01, color = "red", linetype = "dashed", alpha = 0.6, linewidth = 1)+
  labs(x = "Number of mesh cells",
       y = expression("Absolute change in drag coefficient (" * abs(Delta * C[D]) * ")")) +
  pub_theme

# Tag penalty comparison ####

lifecols = c("#EEE7C3", "#86CAC1", "#48748D", "#27305E")

ggplot(cfd_wide_1)+
  geom_line(aes(Shark_size, rel_drag_pct), color = "#EDDA2C", linewidth = 1.2, alpha = 0.8)+
  geom_point(aes(Shark_size, rel_drag_pct), fill = "#EDDA2C", size = 4, shape = 21, alpha = 0.8)+
  geom_line(aes(Shark_size, Tag.Shark_weight), color = "#86CAC1", linewidth = 1.2, alpha = 0.8)+
  geom_point(aes(Shark_size, Tag.Shark_weight), fill = "#86CAC1", size = 4, shape = 21, alpha = 0.8)+
  geom_line(aes(Shark_size, Tag.shark_area), color = "#27305E", linewidth = 1.2, alpha = 0.8)+
  geom_point(aes(Shark_size, Tag.shark_area), fill = "#27305E", size = 4, shape = 21, alpha = 0.8)+
  labs(x = "Shark size (m TL)",
       y = "Tag penalty (%)") +
  pub_theme

# Size ####

ggplot(cfd_wide %>%
         filter(Velocity == 1))+
  geom_line(aes(Shark_size, TOTAL_FORCE_X_NO), linewidth = 2, alpha = 0.2)+
  geom_line(aes(Shark_size, TOTAL_FORCE_X_YES),linewidth = 2, alpha = 0.5)+
  geom_point(aes(Shark_size, TOTAL_FORCE_X_NO), size = 4, alpha = 0.5)+
  geom_point(aes(Shark_size, TOTAL_FORCE_X_YES),size = 4)+
  labs(x = "Shark size (cm TL)",
       y = "Drag Force (N)") +
  pub_theme


# Speed ####
ggplot(cfd_wide %>%
         filter(Shark_size == 1.2))+
  geom_line(aes(Velocity, TOTAL_FORCE_X_NO), linewidth = 1, alpha = 0.2)+
  geom_line(aes(Velocity, TOTAL_FORCE_X_YES),linewidth = 1, alpha = 0.5)+
  geom_point(aes(Velocity, TOTAL_FORCE_X_NO), size = 2, alpha = 0.5)+
  geom_point(aes(Velocity, TOTAL_FORCE_X_YES),size = 2)+
  labs(x = "Swimming speed (m/s)",
       y = "Drag Force (N)") +
  pub_theme
