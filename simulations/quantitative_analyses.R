# Introduction ----
# Title: Role of erythrocyte glycolytic and redox metabolism in exercise
#
# Description: Quantitative analyses using experimental data and mathematical 
# equations from Erythrocyte Metabolism.
#
# Author: Panagiotis N. Chatzinikolaou
# Affiliation: PhD candidate at Aristotle University of Thessaloniki, Greece
# Contact: chatzinpn@phed-sr.auth.gr
#
# Latest update: 13/02/2024

# Functions ----
source(file = "D:/Rdirectory/phd_projects/Erythrocyte-Metabolism/code/RBC_oxygen_model.R")

# Data ----
# Set the experimental parameters

# Experimental 2,3-BPG, pH, hemoglobin and hematocrit values
# Control condition
bpg23_con_rest <- 5.05 * (10 ^ (-3)) 
bpg23_con_pre <- 5.11 * (10 ^ (-3))
bpg23_con_post <- 5.81 * (10 ^ (-3))

hematocrit_con <- 0.429  # 42.9 %
hemoglobin_con <- 149.15 # g/L

blood_pH_con_pre <- 7.42
blood_pH_con_post <- 7.18

# Oxidative stress condition
bpg23_exp_rest <- 5.08 * (10 ^ (-3))
bpg23_exp_pre <- 5.15 * (10 ^ (-3))
bpg23_exp_post <- 6.28 * (10 ^ (-3))

hematocrit_exp <- 0.434  # 43.4 %
hemoglobin_exp <- 149.95 # g/L

blood_pH_exp_pre <- 7.43
blood_pH_exp_post <- 7.17

# Values drawn from the literature (Sun et al 2000)
# We assumed that temperature and pCO2 would be similar between the two conditions

blood_temp_pre <- 37.0 # degrees Celcius
blood_temp_post <- 38.5

blood_pco2 <- 37.9
blood_po2 <- 98.9

# Quantitative analyses ----
# Step 1. Calculate p50 values for each condition and timepoint (pre-, post-exercise).
rbc_pH_con_pre <- rbc_ph(ph = blood_pH_con_pre)
rbc_pH_con_post <- rbc_ph(ph = blood_pH_con_post)
rbc_pH_exp_pre <- rbc_ph(ph = blood_pH_exp_pre)
rbc_pH_exp_post <- rbc_ph(ph = blood_pH_exp_post)

p50_con_pre <- model_p50(dpg_rbc = bpg23_con_pre, 
                          pH_rbc = rbc_pH_con_pre,
                          temp = blood_temp_pre,
                          pco2 = blood_pco2)

p50_con_post <- model_p50(dpg_rbc = bpg23_con_post, 
                          pH_rbc = rbc_pH_con_post,
                          temp = blood_temp_post,
                          pco2 = blood_pco2)

p50_exp_pre <- model_p50(dpg_rbc = bpg23_exp_pre, 
                          pH_rbc = rbc_pH_exp_pre,
                          temp = blood_temp_pre,
                          pco2 = blood_pco2)

p50_exp_post <- model_p50(dpg_rbc = bpg23_exp_post,
                          pH_rbc = rbc_pH_exp_post,
                          temp = blood_temp_post,
                          pco2 = blood_pco2)

# Step 2. Calculate the oxygen saturation using the Hill equation.
po2 <- 0:100  # create a vector of values for the partial pressure of oxygen.

SHbO2_con_pre <- model_hill(po2 = po2, p50 = p50_con_pre)
SHbO2_con_post <- model_hill(po2 = po2, p50 = p50_con_post)

SHbO2_exp_pre <- model_hill(po2 = po2, p50 = p50_exp_pre)
SHbO2_exp_post <- model_hill(po2 = po2, p50 = p50_exp_post)

# Step 3. Plot
library(ggplot2)
# For the plots, we will use the same theme and color palette as in Figures.qmd doc.

# Control condition: pre-post plot
odc_time_con <- 
  ggplot(data = data.frame(po2, SHbO2_con_pre, SHbO2_con_post), 
                       aes(x = po2, y = SHbO2_con_pre*100)) +
  geom_line(color = "black") +
  geom_line(aes(y = SHbO2_con_post*100), color = "black", linetype = "dashed") +
  geom_vline(xintercept = p50_con_pre, linetype = "solid", color = "black") +
  geom_vline(xintercept = p50_con_post, linetype = "dashed", color = "black") +
  geom_segment(aes(x = p50_con_pre, xend = p50_con_post, y = 80, yend = 80),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.3) +
  labs(title = "Control condition: pre- and post-exercise",
       x = "Partial pressure of oxygen (mmHg)",
       y = "Oxygen saturation (%)") + 
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(from = 0, to = 100, by = 20)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(from = 0, to = 100, by = 20)) +
  theme_rbc3() +
  theme(text = element_text(family = "Segoe UI"),
        plot.title = element_text(hjust = 0, size = 14))

# Oxidative stress condition: pre-post plot
odc_time_exp <-
  ggplot(data = data.frame(po2, SHbO2_exp_pre, SHbO2_exp_post), 
         aes(x = po2, y = SHbO2_exp_pre*100)) +
  geom_line(color = "red") +
  geom_line(aes(y = SHbO2_exp_post*100), color = "red", linetype = "dashed") +
  geom_vline(xintercept = p50_exp_pre, linetype = "solid", color = "red") +
  geom_vline(xintercept = p50_exp_post, linetype = "dashed", color = "red") +
  geom_segment(aes(x = p50_exp_pre, xend = p50_exp_post, y = 80, yend = 80),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.3) +
  labs(title = "Oxidative stress condition: pre- and post-exercise",
       x = "Partial pressure of oxygen (mmHg)",
       y = "Oxygen saturation (%)") + 
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(from = 0, to = 100, by = 20)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(from = 0, to = 100, by = 20)) +
  theme_rbc3() +
  theme(text = element_text(family = "Segoe UI"),
        plot.title = element_text(hjust = 0, size = 14))

# Control vs. oxidative stress condition: pre-exercise plot
odc_cond_pre <- 
  ggplot(data = data.frame(po2, SHbO2_con_pre, SHbO2_exp_pre), 
         aes(x = po2, y = SHbO2_con_pre*100)) +
  geom_line(color = "black") +
  geom_line(aes(y = SHbO2_exp_pre*100), color = "red") +
  geom_vline(xintercept = p50_con_pre, linetype = "solid", color = "black") +
  geom_vline(xintercept = p50_exp_pre, linetype = "solid", color = "red") +
#  geom_segment(aes(x = p50_con_pre, xend = p50_exp_pre, y = 80, yend = 80),
#               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.4) +
  labs(title = "Pre-exercise: control and oxidative stress condition",
       x = "Partial pressure of oxygen (mmHg)",
       y = "Oxygen saturation (%)") + 
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(from = 0, to = 100, by = 20)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(from = 0, to = 100, by = 20)) +
  theme_rbc3() +
  theme(text = element_text(family = "Segoe UI"),
        plot.title = element_text(hjust = 0, size = 14))

# Control vs. oxidative stress condition: post-exercise plot
odc_cond_post <- 
  ggplot(data = data.frame(po2, SHbO2_con_post, SHbO2_exp_post), 
         aes(x = po2, y = SHbO2_con_post*100)) +
  geom_line(color = "black", linetype = "dashed") +
  geom_line(aes(y = SHbO2_exp_post*100), color = "red", linetype = "dashed") +
  geom_vline(xintercept = p50_con_post, linetype = "dashed", color = "black") +
  geom_vline(xintercept = p50_exp_post, linetype = "dashed", color = "red") +
  geom_segment(aes(x = p50_exp_post-9.4, xend = p50_exp_post, y = 80, yend = 80),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.3) +
  labs(title = "Post-exercise: control and oxidative stress condition",
       x = "Partial pressure of oxygen (mmHg)",
       y = "Oxygen saturation (%)") + 
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(from = 0, to = 100, by = 20)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(from = 0, to = 100, by = 20)) +
  theme_rbc3() +
  theme(text = element_text(family = "Segoe UI"),
        plot.title = element_text(hjust = 0, size = 14))

# Combine the plots into a panel using patchwork
library(patchwork)

odc_panel <- 
  (odc_time_con + odc_time_exp) /
  (odc_cond_pre + odc_cond_post) + # add annotations
  plot_annotation(tag_levels = "A", tag_suffix = ".")
odc_panel

# Step 4. Export figures
ggsave(plot = odc_panel, filename = "odc_panel.svg", 
       path = "D:/PhD/Papers/Erythrocyte study/Figures/Study design/", 
       width = 26, height = 18, units = "cm", device = "svg", dpi = 320, scale = 1.2)
