# This code is a helper function to calculate summary and inference statistics.

## Define function to calculate summary statistics
## Description:
## The function gives measures of location, dispersion and inference statistics.
## Specifically, it gives the count, mean, standard deviation, standard error of the mean, and 
## confidence interval (default 95%).
##
## Arguments:
##   data: a data frame.
##   measureVar: the name of a column that contains the variable to be summarised
##   ...: use the three dots to add the names of columns that contain grouping variables
##   removeNA: a boolean that indicates whether to ignore NA's
##   alphaLevel: level of statistical significance to calculate CI (default is 95%)
##
##   Adapted from Raincloud plots tutorial github repository: 
##   https://github.com/RainCloudPlots/RainCloudPlots/blob/master/tutorial_R/summarySE.R

summary_stats <- 
  function(data, measureVar, ..., alphaLevel = 0.05, removeNA = FALSE) {
  
    sum_stats <- data |>
      dplyr::group_by(...) |> 
      dplyr::summarise(n = length({{ measureVar }}), 
                       Mean = mean({{ measureVar }}),
                       SD = sd({{ measureVar }}),
                       Median = mean({{ measureVar }}),
                       IQR = IQR({{ measureVar }}),
                       SE = sd({{ measureVar }}) / sqrt(n),
                       Lower_CI = Mean - (qt((1 - 0.05)/2 + 0.5, n - 1) * SE),
                       Upper_CI = Mean + (qt((1 - 0.05)/2 + 0.5, n - 1) * SE),
                       .groups = 'drop')
    
    return(sum_stats)
  }

# Functions to calculate inference statistics
## Define function to calculate the Standard error of the mean.
sem <- function(x) {
  standard_error <- sd(x) / sqrt(length(x))
  return(standard_error)
}

## Define function to calculate the ± 1.96 SEM to the sample mean to construct 
## the upper and lower bounds of the 95% CI.
conf_int_95 <- function(x, n) {
  lowerCI <- mean(x) - (qt((1 - 0.05)/2 + 0.5, n - 1) * sem(x))
  upperCI <- mean(x) + (qt((1 - 0.05)/2 + 0.5, n - 1) * sem(x))
  x_95CI <- c(lowerCI, upperCI)
  return(x_95CI)
}

# Calculate the mean difference and the change from baseline.
mean_diff <- function(x, y) {
  mean_1 <- mean(x)
  mean_2 <- mean(y)
  diff <- mean_2 - mean_1
  return(diff)
}

mean_change <- function(x, y) {
  mean_1 <- mean(x)
  mean_2 <- mean(y)
  mean_diff <- mean_2 - mean_1
  change <- (mean_diff / mean_1) * 100
  return(change)
}