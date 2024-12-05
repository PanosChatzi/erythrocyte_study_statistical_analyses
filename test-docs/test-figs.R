# Introduction ----
# Title: The role of erythrocyte metabolism in performance
# Purpose: R code to create figures
#
# Panagiotis N. Chatzinikolaou
# chatzinpn[at]phed-sr.auth.gr
#
# First created: 03/08/2023
# Last update: 04/10/2023

# Themes ----
# Add default settings

# create a custom palette with your favourite colours. The number of colours must
# match the levels of the factor (e.g., eccentric and concentric condition get two colours).
rbc.palette1 <- c("black", "red")

# add Segoe UI font
windowsFonts(`Segoe UI` = windowsFont('Segoe UI'))

theme_rbc1 <- function() {                                                # creating a new theme function
  ggplot2::theme_bw() +                                                   # using a predefined theme as a base
    theme(axis.text.x = element_text(size = 12, colour = "black"),        # customizing things
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14),
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks.y = element_line(colour = "black"),                  # Change y axis ticks only
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0.5),
          plot.background = element_blank(),
          legend.position = "none",                     # Delete legend
          panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())               # to delete the border use: element_blank()
}                                                       # element_rect(colour = "black", size = 1, fill = NA)

theme_rbc2 <- function() {                                               # creating a new theme function
  ggpubr::theme_classic2() +                                             # using a predefined theme as a base
    theme(axis.text.x = element_text(size = 12, colour = "black"),       # customizing things
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14),
          panel.grid = element_blank(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0.5),
          legend.position = "none") # panel.border = element_blank()
}

theme_rbc3 <- function() {                                               # creating a new theme function
  cowplot::theme_cowplot() +                                             # using a predefined theme as a base
    theme(text = element_text(family = "Segoe UI"),
          legend.position = "none") # panel.border = element_blank()
}

# Packages ----
#source("r_docs/SummaryStats.R")

# Load ggplot2 library
library(tidyr)
library(dplyr)
library(cowplot)
library(ggbeeswarm)
library(ggplot2)

# Read data ----
# Read excel file from start
full_data <- read.csv2("data/master_database_v4.csv")

# Load the saved Rdata files
load("data/data.RData")

# Transform ----
# Convert to tidy format
lac_data <- 
  full_data |> select(ID, lac_con_rest:lac_ecc_post30) |> 
  pivot_longer(!ID,
               names_to = "Time",
               values_to = "lactate") |> 
  separate(Time, into = c("Variable_name", "Condition", "Time"))

# Factors
lac_data$Condition <- factor(lac_data$Condition)
lac_data$ID <- factor(lac_data$ID)

lac_data$Time <- factor(lac_data$Time,
                        levels = c("rest", "pre", "post", "post10", "post30"),
                        labels = c("Rest", "Pre", "Post", "Post-10'", "Post-30'"))

# Summary statistics
(sum_lac <- lac_data |> 
  group_by(Condition, Time) |> 
  summarise(n = length(lactate),
            Mean = mean(lactate),
            Median = mean(lactate),
            SD = sd(lactate),
            SE = sd(lactate) / sqrt(n),
            IQR = IQR(lactate),
            error_margin = qt((1 - 0.05)/2 + 0.5, n - 1) * SE))

# Beeswarm plot ----
# link the custom colours (black and red) with the levels of condition (control, eccentric)
names(rbc.palette1) <- levels(lac_data$Condition)

pd.bee = position_dodge(0.75)
# Dodge.width gives some space between the two conditions, otherwise the points would be 
# on top of each other. See here for more info:
# https://statisticsglobe.com/jitter-position-dodge-simultaneously-ggplot2-plot-r

(plot3 <- 
ggplot(data = lac_data, 
       aes(x = Time, y = lactate, group = Condition, colour = Condition)) + 
  # add raw data points with the ggbeeswarm::quasirandom function
  ggbeeswarm::geom_quasirandom(dodge.width = 0.75,
                               size = 2, width = .15, shape = 1) + #stroke = 0.70
  ggbeeswarm::geom_quasirandom(dodge.width = 0.75, 
                               size = 2, width = .15, alpha = 0.3) +
  # add lines connecting cell means by condition
  stat_summary(fun = "mean", geom = "line", linewidth = 0.7, position = pd.bee) +
  # add 95% CI to cell means
  stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.10, size = 0.70,
                 position = pd.bee) +
  # add data points representing the means across each condition
  stat_summary(fun = "mean", geom = "point", size = 3, position = pd.bee,
               shape = 21, fill = "white") +
  # add custom theme
  theme_rbc3() +
  # add custom labels for 'x' and 'y' axis
  labs(x = "", y = "Glycolysis flux (μM/min)\n") +
  # remove tickmarks from 'x' axis
  ggpubr::rremove("x.ticks"))

# change labels of x axis levels and add custom colours

(plot3 <- plot3 + 
  scale_x_discrete(labels=c("Rest" = "Baseline", 
                            "Pre" = "Pre-arm",
                            "Post" = "Post-arm",
                            "Post-10'" = "Post-10'",
                            "Post-30'" = "Post-30'")) +
  scale_color_manual(values = rbc.palette1))

# PsyTR plot ----
## Nordmann PsyTeachR tutorial
jd <- 0.75
jd2 <- 0.4

(plot4 <- 
  ggplot(data = lac_data,
         aes(x = Time, y = lactate, group = Condition, colour = Condition)) +
  # add raw data points and use dodge.width discern the within conditions
  geom_point(shape = 16, size = 2, alpha = 0.30,
             position = position_jitterdodge(jitter.width = jd2,
                                             dodge.width = jd,
                                             seed = 0)) +
  geom_point(shape = 1, size = 2, #stroke = 0.70
             position = position_jitterdodge(jitter.width = jd2,
                                             dodge.width = jd,
                                             seed = 0)) +  
  # add lines to connect each participant data points across conditions
#  geom_line(aes(group = ID), alpha = 0.2) +
  # add lines connecting cell means by condition
  stat_summary(fun = "mean", geom = "line", linewidth = 0.70,
               position = position_dodge(width = jd)) +
  # add 95% CI to cell means
  stat_summary(fun.data = "mean_cl_normal", position = position_dodge(width = jd),
               geom = "errorbar", width = 0.10, size = 0.70) +
  # add data points representing the means across each condition
  stat_summary(fun = "mean", geom = "point", size = 3, shape = 21, fill = "white",
                position = position_dodge(width = jd)) +
  # add custom theme
  theme_rbc3() +
  # add 'x' and 'y' axis labels
  labs(x = "", y = "Glycolysis flux (μM/min)\n") +
  # remove tick marks in 'x' axis
  ggpubr::rremove("x.ticks"))

plot4 +
  # change names of 'x' axis tick labels
  scale_x_discrete(labels=c("Rest" = "Baseline", 
                            "Pre" = "Pre-arm",
                            "Post" = "Post-arm",
                            "Post-10'" = "Post-10'",
                            "Post-30'" = "Post-30'")) +
  # add manual colour
  scale_color_manual(values = rbc.palette1)

# Facet line ----
facet_names <- c(con = "Control", ecc = "Oxidative stress")

(plot_facet <- 
  ggplot(data = lac_data,
       aes(x = Time, y = lactate, group = ID, colour = Condition)) +
  geom_line(alpha = 0.30) +
  geom_point() +
  facet_wrap(~Condition, labeller = labeller(Condition = facet_names)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "", y = "Glycolysis flux (μM/min)\n") +
  scale_color_manual(values = rbc.palette1) +
  scale_x_discrete(labels=c("Rest" = "Baseline", 
                              "Pre" = "Pre-arm",
                              "Post" = "Post-arm",
                              "Post-10'" = "Post-10'",
                              "Post-30'" = "Post-30'")))

# Facet interaction ----
# create the interaction factors (Time*Condition) and (ID*Condition)
lac_interaction <- lac_data %>%
  mutate(interact = paste(Time, Condition, sep = ", "))

lac_interaction <- lac_interaction |> 
  mutate(intersubject = paste(ID, Condition, sep = ", "))

# sort the interaction levels into the correct order
lac_interaction$interact <- 
  factor(lac_interaction$interact,
         levels = c("Rest, con", "Pre, con", "Post, con", "Post-10', con", "Post-30', con",
                    "Rest, ecc", "Pre, ecc", "Post, ecc", "Post-10', ecc", "Post-30', ecc"))

# create the ggplot with the default aesthetics
(interplot <- 
ggplot(data = lac_interaction,
       aes(x = interact, y = lactate, colour = Condition)) +
  # add individual lines connecting subjects within each condition
  geom_line(alpha = 0.30,
            aes(group = intersubject)) +
  # add individual data points
  geom_point(alpha = 0.70) +
  theme_rbc3() +
  theme(legend.position = "none"))

(interplot <- interplot +
  annotate(geom="text", x = 1, y = 90, label="Control") +
  annotate(geom="text", x = 6.5, y = 90, label="Oxidative stress") +
  scale_x_discrete(labels=c("Baseline", "Pre-arm", "Post-arm", "Post-10'", "Post-30'",
                            "Baseline", "Pre-arm", "Post-arm", "Post-10'", "Post-30'")) +
  scale_color_manual(values = rbc.palette1) +
  # add 'x' and 'y' axis labels
  labs(x = "", y = "Glycolysis flux (μM/min)\n"))

# more concise way for this plot
ggplot(data = lac_data,
       aes(x = interaction(Time, Condition), y = lactate, group = interaction(ID, Condition), colour = Condition)) +
  geom_line(alpha = 0.30) +
  geom_point(alph = 0.7) +
  theme(legend.position = "none") 

# Mean & CI plot ----
(interplot2 <- 
    ggplot(data = lac_data, 
           aes(x = Time, y = lactate, group = Condition, colour = Condition)) + 
    # add lines connecting cell means by condition
    stat_summary(fun = "mean", geom = "line", linewidth = 0.7, position = position_dodge(0.4)) +
    # add 95% CI to cell means
    stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.10, linewidth = 0.70,
                 position = position_dodge(0.4)) +
    # add data points representing the means across each condition
    stat_summary(fun = "mean", geom = "point", size = 3.5, position = position_dodge(0.4),
                 shape = 21, stroke = 0.90, fill = "white") +
    # add custom theme
    theme_rbc3() +
    # add custom labels for 'x' and 'y' axis
    labs(x = "", y = "Glycolysis flux (μM/min)\n") +
    # change labels of x axis levels and add custom colours
    scale_x_discrete(labels=c("Rest" = "Baseline", 
                            "Pre" = "Pre-arm",
                            "Post" = "Post-arm",
                            "Post-10'" = "Post-10'",
                            "Post-30'" = "Post-30'")) +
    scale_color_manual(values = rbc.palette1))

# Combine plots
plot_grid(interplot2, interplot, labels = "auto", ncol = 1, rel_widths = c(1, 2))

# Other plots
# Boxplot - Statistics globe
ggplot(data = lac_data,
       aes(x = Time, y = lactate, fill = Condition)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = 0.75),
             aes(fill = Condition), pch = 21) +
  theme_bw() 

# Mean & CI from R in Action
pd = position_dodge(0.75)
pd2 = position_dodge(0.4)

(plot1 <- 
    ggplot(data = sum_lac, 
           aes(x = Time, y = Mean, group = Condition, color = Condition)) +
    geom_line(linewidth = 1, position = pd2) +
    geom_point(size = 3, position = pd2) +
    geom_errorbar(aes(ymin = Mean - CI, ymax = Mean + CI, width = 0.2), position = pd2, size = 1) +
    #  geom_point(data = lac_data, aes(x = Time, y = lactate, color = Condition), position = pd) +
    theme_rbc1() +
    labs(x = "", y = "Glycolysis flux (μM/min)\n") +
    ggpubr::rremove("x.ticks"))

plot1 + 
  scale_x_discrete(labels=c("Rest" = "Baseline", 
                            "Pre" = "Pre-arm",
                            "Post" = "Post-arm",
                            "Post-10'" = "Post-10'",
                            "Post-30'" = "Post-30'")) +
  scale_color_manual(values = c("red", "black")) 

# Dotplot and summary stats
(plot2 <- 
    ggplot() +
    geom_point(data = lac_data, aes(x = Time, y = lactate, colour = Condition), 
               position = position_jitterdodge(), alpha = 0.4) +
    stat_summary(data = lac_data, aes(x = Time, y = lactate, colour = Condition, group = Condition),
                 fun = mean, geom = "line", linewidth = 1, position = pd) +
    stat_summary(data = lac_data, aes(x = Time, y = lactate, colour = Condition, group = Condition),
                 fun = mean, geom = "point", size = 3, position = pd) +
    geom_errorbar(data = sum_lac, aes(x = Time, ymin = Mean - CI, ymax = Mean + CI, color = Condition),
                  position = pd, width = 0.2, size = 1) +
    theme_rbc1() +
    labs(x = "", y = "Glycolysis flux (μM/min)\n") +
    ggpubr::rremove("x.ticks"))

plot2 +
  scale_x_discrete(labels=c("Rest" = "Baseline", 
                            "Pre" = "Pre-arm",
                            "Post" = "Post-arm",
                            "Post-10'" = "Post-10'",
                            "Post-30'" = "Post-30'")) +
  scale_color_manual(values = c("red", "black")) 

#geom_errorbar(data = sum_lac, aes(x = Time, ymin = Mean - CI, ymax = Mean + CI, color = Condition),
#                position = pd.bee, width = 0.2, linewidth = 1) 

# Save ----
save(lac_data, sum_lac, file = "data/data.RData")

# Save images
# Using ggsave() function
ggsave(plot = glycoplot, filename = "Glycolysis.svg", path = "D:/PhD/Papers/Erythrocyte study/Figures",
       device = "svg", width = 19.82, height = 13.2, units = "cm", dpi = 320)

# Annotations ----

p1 <- ggplot() + geom_point(aes(x=1:100, y=1:100))
p2 <- ggplot() + geom_point(aes(x=1:100, y=1:100))
p3 <- ggplot() + geom_point(aes(x=1:100, y=1:100))
p4 <- ggplot() + geom_point(aes(x=1:100, y=1:100))

(p1 | p2) / (p3 | p4) + plot_annotation(tag_levels = 'A', tag_suffix = '.')

ggarrange(p1, p1, p1, widths = c(3, 1), heights = c(5, 1))

# Create a ggplot and add two titles on the x axis
library(ggplot2)

# Create sample data
data <- data.frame(
  x = c("A", "B", "C", "D"),
  y = c(10, 15, 8, 20)
)

# Create a plot with custom x-axis labels
ggplot(data, aes(x = x, y = y)) +
  geom_point() + labs(x = "") +
  
annotate("text", x = c(2,4.5), y = c(-1, 0), 
 label = c("label 1", "label 2") , color="orange", 
 size=7 , fontface="bold") + coord_cartesian(ylim = c(5, 20), clip = "off")



