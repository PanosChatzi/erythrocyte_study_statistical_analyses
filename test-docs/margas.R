# Packages ----
# Load packages
library(patchwork)
library(tidyr)
library(dplyr)
library(ggbeeswarm)
library(ggplot2)

# Functions ----
# Create your own colour palette (black and red for control and oxidative stress)
rbc.palette1 <- c("black", "red")

# add Segoe UI font family (or your favourite font)
windowsFonts(`Segoe UI` = windowsFont('Segoe UI'))

# Define a custome theme for the plots
## there are numerous available ggplot2 themes, but I use my own
theme_rbc3 <- function() {                                                # creating a new theme function
  theme_bw() + # using a predefined theme as a base
    theme(text = element_text(family = "Segoe UI"),
          axis.text.x = element_text(size = 12, colour = "black"), 
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14),
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(colour = "black"), # Change y axis ticks only
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0.5),
          plot.background = element_blank(),
          legend.position = "none",    # Delete legend
          panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())  # to delete the border use: element_blank()
} # element_rect(colour = "black", size = 1, fill = NA)

# Data ----
# Load data
# For RData fles we use load() function
load("data/tidyData.RData")

# If you have a CSV file, you can use read.csv() (base R) or read_csv() from readr package

# RBC figures ----
# set custom dodge
bee_width = 0.15 # spread the bees (width of the data points)
bee_dodge = 0.75 # dodge the bees (space between the groups)

# Create a plot for the lactate flux
(glycoplot <- 
    lactate_dat_tidy |> 
    ggplot(aes(x = timepoint, y = glycolysis_flux, group = condition, colour = condition)) + 
    # add raw data points with the ggbeeswarm::quasirandom function
    ggbeeswarm::geom_quasirandom(dodge.width = bee_dodge, width = bee_width,
                                 size = 2.5, shape = 1) + #stroke = 0.70
    ggbeeswarm::geom_quasirandom(dodge.width = bee_dodge, width = bee_width,
                                 size = 2.5, alpha = 0.3) +
    # add lines connecting cell means by condition
    stat_summary(fun = "mean", geom = "line", linewidth = 0.70, 
                 position = position_dodge(bee_dodge)) +
    # add 95% CI to cell means
    stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.10, 
                 linewidth = 0.70, position = position_dodge(bee_dodge)) +
    # add data points representing the means across each condition
    stat_summary(fun = "mean", geom = "point", size = 3.2, stroke = 1,
                 position = position_dodge(bee_dodge), shape = 21, fill = "white") +
    # add custom theme
    theme_rbc3() +
    # add custom labels for 'x' and 'y' axis
    labs(x = "", y = "Glycolytic flux (μM/min)") +
    # change 'x' axis labels
    scale_x_discrete(labels=c("Baseline" = "Baseline", 
                              "Pre-arm" = "Pre-exercise",
                              "Post-arm" = "Post-exercise",
                              "Post-10'" = "Post-10' exercise",
                              "Post-30'" = "Post-30' exercise")) +
    # add custom colours
    scale_color_manual(values = rbc.palette1))

# Change axes limits
(glycoplot <- glycoplot + 
    coord_cartesian(ylim = c(0, 100)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(from = 0, to = 100, by = 20)))

# Save
ggsave(plot = glycoplot, # choose the plot you want to save
       filename = "Glycolysis.svg", # choose file name
       path = "Figures/",   # choose folder path
       device = "svg",      # figure format (png, svg, pdf, etc)
       width = 22.31, height = 13.42, units = "cm", dpi = 320) # figure size and resolution


# Try again for vo2peak
# Arms
vo2_arm_plot <- 
  vo2_dat_tidy |> 
  ggplot(aes(x = condition, y = vo2Arm, group = 1, colour = condition)) + 
  ggbeeswarm::geom_quasirandom(width = 0.06, dodge.width = NULL,
                               size = 2.5, shape = 1) + #stroke = 0.70
  ggbeeswarm::geom_quasirandom(width = 0.06, dodge.width = NULL,
                               size = 2.5, alpha = 0.3) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.70,
               colour = "black") + # add a separate colour
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.04,
               linewidth = 0.70, colour = "black") +
  stat_summary(fun = "mean", geom = "point", size = 3.2, stroke = 1,
               shape = 21, fill = "white", colour = "black") +
  theme_rbc3() +
  labs(x = "", y = expression("Arms VO"[2]~"peak (ml/min/kg)")) +
  scale_color_manual(values = rbc.palette1)

(vo2_arm_plot <- 
    vo2_arm_plot + 
    coord_cartesian(ylim = c(0, 60)) + 
    scale_y_continuous(expand = c(0, 0), 
                       breaks = seq(from = 0, to = 60, by = 10)))

# Legs
vo2_leg_plot <- 
  vo2_dat_tidy |> 
  ggplot(aes(x = condition, y = vo2Leg, group = 1, colour = condition)) + 
  ggbeeswarm::geom_quasirandom(width = 0.06, dodge.width = NULL,
                               size = 2.5, shape = 1) + #stroke = 0.70
  ggbeeswarm::geom_quasirandom(width = 0.06, dodge.width = NULL,
                               size = 2.5, alpha = 0.3) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.70,
               colour = "black") +
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.04,
               linewidth = 0.70, colour = "black") +
  stat_summary(fun = "mean", geom = "point", size = 3.2, stroke = 1,
               shape = 21, fill = "white", colour = "black") +
  theme_rbc3() +
  labs(x = "", y = expression("Legs VO"[2]~"peak (ml/min/kg)")) +
  scale_color_manual(values = rbc.palette1)

(vo2_leg_plot <- 
    vo2_leg_plot + 
    coord_cartesian(ylim = c(0, 60)) + 
    scale_y_continuous(expand = c(0, 0), 
                       breaks = seq(from = 0, to = 60, by = 10)))

# Combine and save the plots (check patchwork for guides)
(vo2_panel <- vo2_arm_plot + vo2_leg_plot + 
    plot_annotation(tag_levels = 'A', tag_suffix = '.'))

ggsave(plot = vo2_panel, filename = "vo2_panel.svg", 
       path = "D:/PhD/Papers/Erythrocyte study/Figures/Physiology/", 
       width = 26.43, height = 11.09, units = "cm", device = "svg", dpi = 320)

