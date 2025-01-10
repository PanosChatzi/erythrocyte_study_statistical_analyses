# Introduction ----
# Title: Role of erythrocyte glycolytic and redox metabolism in exercise
#
# Description: Calculate p50 values of each participant at each condition and timepoint.
#
# Author: Panagiotis N. Chatzinikolaou
# Affiliation: PhD candidate at Aristotle University of Thessaloniki, Greece
# Contact: chatzinpn@phed-sr.auth.gr
#
# Latest update: 11/09/2024

# Functions ----
source(file = "D:/Rdirectory/phd_projects/Erythrocyte-Metabolism/code/RBC_oxygen_model.R")

# Data ----
p50_data <- readr::read_csv2("data/p50_data.csv")

blood_mean_df <- colMeans(p50_data[, 2:ncol(p50_data)])
blood_mean_df

# Temperature values were drawn from the literature (Sun et al 2000)
# We assumed that temperature, pO2 pCO2 would be similar between the two conditions
blood_temp_pre <- 37.0 # degrees Celcius
blood_temp_post <- 38.5

blood_pco2 <- 40
blood_po2 <- 100

# Modeling ----
# Step 1. Calculate p50 values for each condition and timepoint (pre-, post-exercise).
# Control condition p50
p50_con_rest <- model_p50(dpg_rbc = p50_data$bpg_con_rest,
                          pH_rbc = rbc_ph(p50_data$pH_con_rest),
                          temp = blood_temp_pre,
                          pco2 = blood_pco2)

p50_con_pre <- model_p50(dpg_rbc = p50_data$bpg_con_pre, 
                         pH_rbc = rbc_ph(p50_data$pH_con_pre),
                         temp = blood_temp_pre,
                         pco2 = blood_pco2)

p50_con_post <- model_p50(dpg_rbc = p50_data$bpg_con_post, 
                          pH_rbc = rbc_ph(p50_data$pH_con_post),
                          temp = blood_temp_post,
                          pco2 = blood_pco2)

# Experimental condition p50
p50_exp_rest <- model_p50(dpg_rbc = p50_data$bpg_exp_rest,
                          pH_rbc = rbc_ph(p50_data$pH_exp_rest),
                          temp = blood_temp_pre,
                          pco2 = blood_pco2)

p50_exp_pre <- model_p50(dpg_rbc = p50_data$bpg_exp_pre, 
                         pH_rbc = rbc_ph(p50_data$pH_exp_pre),
                         temp = blood_temp_pre,
                         pco2 = blood_pco2)

p50_exp_post <- model_p50(dpg_rbc = p50_data$bpg_exp_post,
                          pH_rbc = rbc_ph(p50_data$pH_exp_post),
                          temp = blood_temp_post,
                          pco2 = blood_pco2)

# Step 2. Prepare dataframe
p50_df <- data.frame(p50_con_rest, p50_con_pre, p50_con_post, 
                     p50_exp_rest, p50_exp_pre, p50_exp_post)

p50_df$ID <- 1:nrow(p50_df)
p50_df <- p50_df[,c(7, 1:6)]

p50_tidy <- 
  p50_df |> 
  tidyr::pivot_longer(cols = -ID,
                      cols_vary = "fastest",
                      names_to = c(".value", "condition", "timepoint"), 
                      names_sep = "_")

# Calculate the mean values
p50_means <- colMeans(p50_df[, 2:ncol(p50_df)])
p50_means

# Step 4. Save data in RData file
save(p50_tidy, file = "quantitative_analysis/quantitative results/p50_model_results.RData")
#write.csv(p50_df, file = "data/p50_model_results.csv")

# Statistics ----
# Statistical analysis for p50 data
library(afex)     # Anova
library(emmeans)  # Estimated marginal means and post-hocs
library(multcomp) # Multiple pairwise comparisons

# Step 1. Fit anova and lm models
# Fit aov with afex
p50_aov_res <- 
  aov_ez(id = "ID",
         dv = "p50",
         data = p50_tidy, 
         within = c('condition', 'timepoint'),
         anova_table = list(es = "pes"))

p50_aov_nice <- knitr::kable(nice(p50_aov_res))
p50_aov_nice

# Fit lm model
# Finally, we fit a linear model for later use (see post-hoc effect sizes section)
lm_model <- lm(p50 ~ condition * timepoint,
               data = p50_tidy,
               contrasts = list(condition = contr.sum, timepoint = contr.sum))

sigma(lm_model)
sd(residuals(p50_aov_res))

# Step 2. Post-hocs with emmeans and multcomp
emms_obj <- emmeans(p50_aov_res, ~ condition * timepoint)

# Use pairs for the pairwise comparisons
pairs_obj <- pairs(emms_obj, adjust = "bonferroni")
post_hocs_obj <- test(pairs_obj)

summary(as.glht(pairs_obj), test = adjusted("bonferroni"))

# Step 3. Effect sizes and clean
# Calculate Hedge's g correction factor J
n <- 20       # Sample size
df <- n - 1   # Degrees of freedom for paired samples
# Using Cousineau & Goulet-Pelletier equation
J <- exp ( lgamma(df / 2) - log(sqrt(df / 2)) - lgamma((df - 1) / 2) )

# Effect sizes calculation
effsize_obj <- eff_size(emms_obj,
                        sigma = (sigma(lm_model)),
                        edf = df.residual(lm_model))

effsize_df <- as.data.frame(effsize_obj,
                            row.names = NULL,
                            check.names = TRUE,
                            destroy.annotations = TRUE)

post_hocs_df <- as.data.frame(post_hocs_obj,
                              row.names = NULL,
                              check.names = TRUE,
                              destroy.annotations = TRUE)

effsize_df$p.values <- post_hocs_df$p.value

# Create a new column with the hedge's g values
effsize_df$hedgesG <- effsize_df$effect.size * J

# Clean
# Subset based on p-values and remove NA values
effsize_df_2 <- na.omit(effsize_df[effsize_df$p.values < 0.05, ])

# Add a new column with asterisks based on significance
effsize_df_2$asteriscs <- 
  ifelse(effsize_df_2$p.values < 0.001, "***", 
         ifelse(effsize_df_2$p.values < 0.01, "**", 
                ifelse(effsize_df_2$p.values < 0.05, "*", "")))

# Keep only the columns of interest for conciseness
effsize_df_3 <- subset(effsize_df_2, 
                       select = c(contrast, hedgesG, p.values, asteriscs))

# Finally, sort by alphabetical order for ease of read
effsize_df_final <- effsize_df_3[order(effsize_df_3$contrast), ]
effsize_df_final

p50_results <- 
  list(Summary_table = p50_means,
       aov_summary = p50_aov_nice, 
       post_hocs = post_hocs_obj, 
       effect_sizes = effsize_obj, 
       sign_es = effsize_df_final)

# Save
save(p50_results, file = "quantitative_analysis/quantitative results/p50_statistics.RData")

# pH ANOVA ----
# Tidy the data
pH_df <- p50_data[, 11:ncol(p50_data)]

pH_df$ID <- 1:nrow(pH_df)

pH_df <- pH_df[,c(7, 1:6)]
pH_df

pH_tidy_dat <- 
  pH_df |> 
  tidyr::pivot_longer(cols = -ID,
                      cols_vary = "fastest",
                      names_to = c(".value", "condition", "timepoint"), 
                      names_sep = "_")
pH_tidy_dat

# Fit aov model
ph_aov_res <- 
aov_ez(id = "ID",
       dv = "pH",
       data = pH_tidy_dat, 
       within = c('condition', 'timepoint'),
       anova_table = list(es = "pes"))

ph_aov_nice <- knitr::kable(nice(ph_aov_res))
ph_aov_nice

#
load(file = "data/nirs_arms_occl_tidy.RData")

arm_thb_aov_res <- 
  afex::aov_ez(id = "ID",
               dv = "tHb",
               data = nirs_arms_occl_tidy, 
               within = c('condition', 'time'),
               anova_table = list(es = "pes"))
knitr::kable(afex::nice(arm_thb_aov_res))
