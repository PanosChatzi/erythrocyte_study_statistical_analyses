# Custom functions for calculating Cohen's d for paired samples
# The custom functions were derived from the Guide to Effect Sizes and Confidence Intervals 
# (2023) and Lakens (2013). Find links at the end of the file.

# Calculate the effect sizes of paired sample(s).
# A. Cohen's Dz: Standardized mean difference (change score)
## A.1. Raw data
sdm_z <- function(x, y, r) {
  
  mean_1 <- mean(x)
  mean_2 <- mean(y)
  mean_diff <- mean_2 - mean_1
  
  sd_1 <- sd(x)
  sd_2 <- sd(y)
  sum_sqr_sd <- (sd_1 ^ 2 + sd_2 ^ 2)
  sd_diff <- sqrt(sum_sqr_sd - 2 * r * sd_1 * sd_2)
  
  sdmdiff <- mean_diff / sd_diff
  
  return(sdmdiff)
}

## A.2. Summary data
sdm_z_sum <- function(mean_1, mean_2, sd_1, sd_2, r) {
  
  mean_diff <- mean_2 - mean_1
  
  sum_sqr_sd <- (sd_1 ^ 2 + sd_2 ^ 2)
  sd_diff <- sqrt(sum_sqr_sd - 2 * r * sd_1 * sd_2)
  
  sdmdiff <- mean_diff / sd_diff
  
  return(sdmdiff)
}

# B. Cohen's Drm: repeated measures data
# 1. Using raw data
sdm_rm <- function(x, y, r = 0.5) {
  
  mean_diff <- mean(y) - mean(x)
  
  sd_1 <- sd(x)
  sd_2 <- sd(y)
  sum_sqr_sd <- (sd_1 ^ 2 + sd_2 ^ 2)
  sd_diff <- sqrt(sum_sqr_sd - 2 * r * sd_1 * sd_2)
  
  sdm <- (mean_diff / sd_diff) * sqrt(2 * (1 - r)) # Based on Lakens (2013).
  
  return(sdm)
}

# 2. Using summary data
sdm_rm_sum <- function(mean_1, mean_2, sd_1, sd_2, r = 0.5) {
  
  mean_diff <- mean_2 - mean_1
  
  sum_sqr_sd <- (sd_1 ^ 2 + sd_2 ^ 2)
  sd_diff <- sqrt(sum_sqr_sd - 2 * r * sd_1 * sd_2)
  
  sdm <- (mean_diff / sd_diff) * sqrt(2 * (1 - r)) # Based on Lakens (2013).
  return(sdm)
}

# C. Cohen's davg: Average variance d
## C.1. Raw data
sdm_avg <- function(x, y, n, CI95 = FALSE) {
  
  mean_1 <- mean(x)
  mean_2 <- mean(y)
  mean_diff <- mean_2 - mean_1
  
  sd_1 <- sd(x)
  sd_2 <- sd(y)
  sd_pooled <- sqrt(((sd_1^2 + sd_2^2) / 2))
  
  davg <- mean_diff / sd_pooled
  
  if(CI95 == TRUE) {
    
    sem_avg <- sqrt( (2 / n) + ((davg^2) / (4 * n)))
    
    ci_low <- davg - (1.96 * sem_avg)
    ci_upp <- davg + (1.96 * sem_avg)
    
    res <- c("davg" = davg, "CI_low" = ci_low, "CI_upp" = ci_upp)
    return(res)
    
  } else {
    
    return(davg)
  }
}

## C.2. Summary data
sdm_avg_sum <- function(mean_1, mean_2, sd_1, sd_2) {
  
  mean_diff <- mean_2 - mean_1
  sd_pooled <- sqrt(((sd_1^2 + sd_2^2) / 2))
  
  davg <- mean_diff / sd_pooled
  return(davg)
}

# D. Becker's db
# This uses the baseline standard deviation when the comparison is a pre/post design
## D.1. Raw data
sdm_b <- function(x, y) {
  
  mean_diff <- mean(y) - mean(x)
  sd_pre <- sd(x)
  
  sdmb <- mean_diff / sd_pre
  return(sdmb)
}

## D.2. Summary data
sdm_b_sum <- function(mean_pre, mean_post, sd_pre) {
  
  mean_diff <- mean_post - mean_pre
  
  sdmb <- mean_diff / sd_pre
  return(sdmb)
}

# E. Convert the t statistic to Cohen's dz
sdm_t <- function(t, n) {
  sdm_t <- t / sqrt(n)
  return((sdm_t))
}

# Standard error of the SDM repeated measures.
se_sdm <- function(sdm, n, r){
  
  se <- sqrt((1 / n) + ((sdm ^ 2) / (2 * n))) * sqrt(2 * (1 - r))
  
  return(se)
}

# References
# Matthiew B. Jane et al. 2024. Guide to Effect Sizes and Confidence Intervals. 
# https://matthewbjane.quarto.pub/guide-to-effect-sizes-and-confidence-intervals

# Lakens. 2013. Calculating and reporting effect sizes to facilitate cumulative science:
# a practical primer for t-tests and ANOVAs.
# doi: 10.3389/fpsyg.2013.00863