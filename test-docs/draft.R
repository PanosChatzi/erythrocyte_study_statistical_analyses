source("r_docs/MyStatsFunctions.R")
source("r_docs/MyCohensEffSizes.R")

# Calculate the summary statistics
sum_dat <- summary_stats(data = performance_tidy, 
                         measureVar = isoTorLeg,
                         condition, timepoint)

mean_pre <- sum_dat$Mean[3]
mean_post <- sum_dat$Mean[4]
sd_pre <- sum_dat$SD[3]
sd_post <- sum_dat$SD[4]

r_1 <- subset(x = performance_tidy,
              subset = condition == "Oxidative stress" & timepoint == "Baseline",
              select = isoTorLeg, 
              drop = TRUE)

r_2 <- subset(x = performance_tidy,
              subset = condition == "Oxidative stress" & timepoint == "Post-48h",
              select = isoTorLeg, 
              drop = TRUE)

corr_sdm <- cor(r_1, r_2)

# Cohen's D avg
sdm_avg_sum(mean_1 = mean_pre,
            mean_2 = mean_post,
            sd_1 = sd_pre,
            sd_2 = sd_post)

sdm_z_sum(mean_1 = mean_pre,
          mean_2 = mean_post,
          sd_1 = sd_pre,
          sd_2 = sd_post,
          r = 0.824)

sdm_rm_sum(mean_1 = mean_pre,
           mean_2 = mean_post,
           sd_1 = sd_pre,
           sd_2 = sd_post,
           r = 0.824)

sdm_b_sum(mean_pre = mean_pre,
          mean_post = mean_post,
          sd_pre = sd_pre)

# Dz with MOTE
sd_diff_mote <- sqrt((sd_pre^2 + sd_post^2) - (2 * corr_sdm * sd_pre * sd_post))

stats <- MOTE::d.dep.t.diff(
  m = mean_post - mean_pre,
  sd = sd_diff_mote,
  n = 20,
  a = 0.05)

# Dz with MOTE
sd_diff_mote <- sqrt((sd_pre^2 + sd_post^2) - (2 * corr_sdm * sd_pre * sd_post))

stats <- MOTE::d.dep.t.diff(
  m = mean_post - mean_pre,
  sd = sd_diff_mote,
  n = 20,
  a = 0.05)

# Drm with MOTE
stats <- MOTE::d.dep.t.rm(
  m1 = mean_pre,
  m2 = mean_post,
  sd1 = sd_pre,
  sd2 = sd_post,
  r = 0.73,
  n = 20,
  a = 0.05
)

# Davg with MOTE
stats <- 
  MOTE::d.dep.t.avg(m1 = mean_pre,
                    m2 = mean_post,
                    sd1 = sd_pre,
                    sd2 = sd_post,
                    n = 20,
                    a = 0.05)

data.frame(d = MOTE::apa(stats$d), 
           dlow = MOTE::apa(stats$dlow), 
           dhigh = MOTE::apa(stats$dhigh))

# Calculate the pooled SD
temp_dat <- metabolic_tidy |> 
  select(ID, hexokinase, condition, timepoint) |> 
  group_by(condition, timepoint) |> 
  summarise(sd_values = sd(hexokinase), .groups = 'drop') |> 
  mutate(sd_names = c("sd1", "sd2", "sd3", "sd4", "sd5", "sd6")) |>
  select(sd_names, sd_values) |> 
  pivot_wider(names_from = sd_names, values_from = sd_values)
# transmute

sqrt(sum(temp_dat^2)/6)


# Helper functions to convert glht and emmGrid objects
# glht objects to tidy objects
# https://rdrr.io/cran/broom/man/tidy.summary.glht.html

# Claas Heuer, September 2015
# https://gist.github.com/cheuerde/3acc1879dc397a1adfb0
glht.table <- function(x) {
  
  pq <- summary(x)$test
  mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
  error <- attr(pq$pvalues, "error")
  pname <- switch(x$alternativ, less = paste("Pr(<", ifelse(x$df ==0, "z", "t"), ")", sep = ""), 
                  greater = paste("Pr(>", ifelse(x$df == 0, "z", "t"), ")", sep = ""), two.sided = paste("Pr(>|",ifelse(x$df == 0, "z", "t"), "|)", sep = ""))
  colnames(mtests) <- c("Estimate", "Std. Error", ifelse(x$df ==0, "z value", "t value"), pname)
  return(mtests)
  
}

emmeans.dataframe <- function(x) {
  
  emm <- summary(x)
  emm_dat <- as.data.frame(emm[])
  return(emm_dat)
  
}
