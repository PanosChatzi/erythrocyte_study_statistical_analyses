# Introduction ----
# Title: Testing statistical analyses
#
# Description: Testing R code for stats
#
# Panagiotis N. Chatzinikolaou
# chatzinpn[at]phed-sr.auth.gr
#
# First created: 15/12/2023
# Last update: 15/12/2023

# Packages ----
source("r_docs/MyStatsFunctions.R")

# Load required libraries
library(dplyr)

# Data ----
# Load the saved Rdata files
load("data/tidyData.RData")

# Stats ----

# A. Glycolysis

# Working examples using Hexokinase
#subset(metabolic_tidy,
#       select = c(ID, condition, timepoint, hexokinase))

# Descriptive statistics
# Calculate the descriptive statistics for the hexokinase variable.
hexokinase_descriptive <- summary_stats(hexokinase_dat, hexokinase, condition, timepoint)

# T-Tests ----
# Shapiro-Wilk test for normality
diff_fun <- function(x, y) {diff = x - y; return(diff)}
arm_vo2_diff <- 
  vo2_dat_tidy[vo2_dat_tidy$condition == "Control", "vo2Arm"] - vo2_dat_tidy[vo2_dat_tidy$condition == "Oxidative stress", "vo2Arm"]

mean_pre <- subset(performance_tidy, 
                   subset = condition == "Oxidative stress" & timepoint == "Baseline", 
                   select = isoTorLeg,
                   drop = TRUE)

mean_post <- subset(performance_tidy, 
                    subset = condition == "Oxidative stress" & timepoint == "Post-48h", 
                    select = isoTorLeg,
                    drop = TRUE)

leg_vo2_diff <- diff_fun(subset(vo2_dat_tidy, condition == "Control", vo2Leg),
                         subset(vo2_dat_tidy, condition == "Oxidative stress", vo2Leg))

shapiro.test(arm_vo2_diff$vo2Arm)

t.test(x = mean_pre,
       y = mean_post,
       alternative = "two.sided",
       mu = 0, 
       paired = TRUE,   
       var.equal = TRUE,
       conf.level = 0.95)

with(vo2_dat_tidy,
  t.test(vo2Arm ~ condition,
       alternative = "two.sided",
       mu = 0, 
       paired = TRUE,   
       var.equal = TRUE,
       conf.level = 0.95))

with(vo2_dat_tidy,
     wilcox.test(vo2Leg ~ condition,
            alternative = "two.sided",
            mu = 0, 
            paired = TRUE,   
            conf.level = 0.95))

effsize::cohen.d(hrLeg ~ condition | Subject(ID),
                 data = vo2_dat_tidy,
                 pooled = TRUE,
                 paired = TRUE,
                 hedges.correction = TRUE)

###########################
mean_pre <- mean(vo2_dat_tidy[vo2_dat_tidy$condition == "Control", "vo2Leg", drop = TRUE])
mean_post <-  mean(vo2_dat_tidy[vo2_dat_tidy$condition == "Oxidative stress", "vo2Leg", drop = TRUE])
sd_pre <- sd(vo2_dat_tidy[vo2_dat_tidy$condition == "Control", "vo2Leg", drop = TRUE])
sd_post <- sd(vo2_dat_tidy[vo2_dat_tidy$condition == "Oxidative stress", "vo2Leg", drop = TRUE])

# Davg with MOTE
stats <- 
  MOTE::d.dep.t.avg(m1 = mean_pre,
                    m2 = mean_post,
                    sd1 = sd_pre,
                    sd2 = sd_post,
                    n = 20,
                    a = 0.05)

# Dz with MOTE
sd_diff_mote <- sqrt((sd_pre^2 + sd_post^2) - (2 * 0.824 * sd_pre * sd_post))

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
  r = 0.824,
  n = 20,
  a = 0.05)

data.frame(d = MOTE::apa(stats$d), 
           dlow = MOTE::apa(stats$dlow), 
           dhigh = MOTE::apa(stats$dhigh))

effsize::cohen.d(vo2Leg ~ condition | Subject(ID),
                 data = vo2_dat_tidy,
                 pooled = TRUE,
                 paired = TRUE,
                 hedges.correction = F)

TOSTER::smd_calc(vo2Leg ~ condition,
                 data = vo2_dat_tidy,
                 paired = TRUE,
                 var.equal = TRUE,
                 bias_correction = F,
                 rm_correction = TRUE)

effectsize::cohens_d(x = vo2_dat_tidy[vo2_dat_tidy$condition == "Control", "vo2Arm", drop = TRUE],
                     y = vo2_dat_tidy[vo2_dat_tidy$condition == "Oxidative stress", "vo2Arm", drop = TRUE],
                     data = vo2_dat_tidy,
                     paired = TRUE,
                     pooled_sd = TRUE)

effectsize::hedges_g(x = vo2_dat_tidy[vo2_dat_tidy$condition == "Control", "vo2Arm", drop = TRUE],
                     y = vo2_dat_tidy[vo2_dat_tidy$condition == "Oxidative stress", "vo2Arm", drop = TRUE],
                     data = vo2_dat_tidy,
                     paired = TRUE,
                     pooled_sd = TRUE)

# Similar results to JASP/Jamovi
rstatix::cohens_d(vo2Leg ~ condition, 
                  data = vo2_dat_tidy,
                  paired = TRUE,
                  var.equal = TRUE,
                  hedges.correction = TRUE)

rstatix::wilcox_effsize(vo2Leg ~ condition,
                        data = vo2_dat_tidy,
                        paired = TRUE,
                        ci = TRUE)

set.seed(2124)
TOSTER::ses_calc(vo2Leg ~ condition,
                 data = vo2_dat_tidy,
                 paired = TRUE)

# Anova ----
## a. Base R ----

# aov() function
# Fit a two-way anova with base R using condition and timepoint as both within factors
hexokinase_aov <- aov(hexokinase ~ condition * timepoint + Error(ID/(condition * timepoint)), 
                      data = metabolic_tidy)

summary(hexokinase_aov)

# Use the lm() function.
options(contrasts = c("contr.helmert", "contr.poly"))

hexokinase_lm <- lm(hexokinase ~ condition * timepoint,
  data = metabolic_tidy)
#  contrasts = list(condition = contr.sum, timepoint = contr.sum))

# This seems to provide the same results as the GAMLj module in Jamovi.

## b. ez ----
# Using the ez package
library(ez)
# Fit a linear mixed model with ez using condition and timepoint as both within factors
# and as their interaction. # Include the partial eta square effect sizes.
hexokinase_ez <- ezANOVA(data = metabolic_tidy, 
                         dv = hexokinase, 
                         wid = ID, 
                         within = c(condition, timepoint), 
                         detailed = TRUE, 
                         type = 3)
hexokinase_ez

# This works! Hurrah!

## c. afex ----
# Using the afex package
library(afex)

# When a factor has more than two levels, call afex options to specify the sphericity correction.
#afex_options(correction_aov = "GG") # Defaut is "GG"

# Fit a linear mixed model with afex using condition and timepoint as both within factors
# and as their interaction. # Include the partial eta square effect sizes.
hexokinase_aov2 <- aov_ez(id = "ID",        # The id argument specifies the subject identifier.
                          dv = "hexokinase", # The dv argument specifies the dependent variable.
                          data = metabolic_tidy, # The data argument specifies the data frame.
                          within = c('condition', 'timepoint'), # The within argument specifies the within-subjects factors.
                          anova_table = list(es = "pes")) # correction = "GG". Specify  effect sizes and sphericity corrections.

hexokinase_aov2
summary(hexokinase_aov2)
# Use the nice function to get the ANOVA table.
nice(hexokinase_aov2)

(model_car <- aov_car(isoTorLeg ~ condition * timepoint + Error(ID/(condition * timepoint)), 
                     data = performance_tidy,
                     anova_table = list(es = 'pes')))

(model_ez <- aov_ez(id = "ID",        # The id argument specifies the subject identifier.
                   dv = "isoTorLeg", # The dv argument specifies the dependent variable.
                   data = performance_tidy, # The data argument specifies the data frame.
                   within = c('condition', 'timepoint'), # The within argument specifies the within-subjects factors.
                   anova_table = list(es = "pes")))# correction = "GG". Specify  effect sizes and sphericity corrections.

# Eta squared ----
# Using the effectsize package
library(effectsize)
# https://github.com/easystats/effectsize

# Partial eta squared
# Use the partial_eta_squared function to get the partial eta squared effect sizes.
eta_squared(hexokinase_aov, ci = 0.90, alternative = "two.sided")     # aov model
eta_squared(hexokinase_aov2, ci = 0.90, alternative = "two.sided")    # afex model

omega_squared(hexokinase_aov, ci = 0.90)                      # aov model
omega_squared(hexokinase_aov2, ci = 0.90)                     # afex model

# Using the MBESS package
library(MBESS)

Lims <- conf.limits.ncf(F.value = 7.3, 
                        conf.level = 0.90, 
                        df.1 <- 1, 
                        df.2 <- 19)

Lower.lim <- Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Upper.lim <- Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)
Lower.lim
Upper.lim

# Using apaTables package
library(apaTables)

get.ci.partial.eta.squared(F.value = 7.3,
                           df1 = 1,
                           df2 = 19,
                           conf.level = 0.90)

# Post-hocs ----
## a. Tukey ----
hexokinase_tukey <- TukeyHSD(aov(hexokinase ~ condition * timepoint, data = hexokinase_dat),
                             which = "condition:timepoint")
## b. emmeans ----
library(emmeans)

emmeans(hexokinase_aov2, pairwise ~ condition * timepoint)
emmeans(hexokinase_aov2, ~ condition|timepoint)

hk_emms <- emmeans(hexokinase_aov2, ~ condition * timepoint)
emmip(hk_emms, condition ~ timepoint, CIs = TRUE)
hk_pairs <- pairs(hk_emms, adjust = "sidak")

confint(pairs(hk_emms, adjust = "sidak"))
test(pairs(hk_emms, adjust = "sidak"))

contrast(hk_emms, "pairwise", adjust = "sidak")
contrast(pairs(hk_emms), 
         interaction = c("pairwise", "consec"), 
         adjust = "bonferroni")

hk_pairs_df <- as.data.frame(hk_pairs,
                             row.names = NULL,
                             check.names = TRUE,
                             destroy.annotations = TRUE)

## Using the aov() model gives the same p values as Jasp.
## Using the afex() model gives the same p values as Jamovi.

## c. multcomp ----
# Use the multcomp package
library(multcomp)
library(multcompView)

# cld provides a compact letter display of the results.
cld(hk_emms, Letters=LETTERS, adjust = "bonferroni")

# Run post-hoc tests using the as.glht function.
summary(as.glht(pairs(emms_obj)), test=adjusted("bonferroni"))

# Display the adjusted results
print(adjusted_posthoc)

# Effect sizes ----
# Find the standard deviation of the models
sigma(hexokinase_lm)
df.residual(hexokinase_lm)

(aov2_sd <- sd(residuals(hexokinase_aov2)))

(sig <- sqrt(hexokinase_aov2$anova_table["condition:timepoint", "MSE"]))
(sig.df <- hexokinase_aov2$anova_table["condition:timepoint", "den Df"])

sigma(hexokinase_lmer)

sd(hexokinase_lm$residuals)

resid_ms <- sum(hexokinase_aov$ID$residuals^2)/hexokinase_aov$ID$df.residual
aov_sd <- sd(hexokinase_aov$ID$residuals)

# Extract the degrees of freedom
df.residual(hexokinase_lmer)
hexokinase_aov$ID$df.residual
hexokinase_aov2$Anova$error.df

# Calculate the post-hoc effect sizes
hk_effsize_lm <- eff_size(hk_emms,
                          sigma = sigma(hexokinase_lm)*sqrt(2),
                          edf = df.residual(hexokinase_lm))

hk_effsize_aov <- 

eff_size(pairs(hk_emms),
         sigma = sigma(hexokinase_lm),
         edf = df.residual(hexokinase_lm),
         method = "identity")

hk_effsize <- as.data.frame(hk_effsize,
                            row.names = NULL,
                            check.names = TRUE,
                            destroy.annotations = TRUE)

hk_effsize |> 
  mutate(hedgesG = effect.size*J)

# Hedge's g correction factor J
n <- 20
df <- n-1
# Using Borenstein & Lakens equation
(J <- (1 - (3 / (4 * df - 1))))

# Using Cousineau & Goulet-Pelletier equation
(J <- exp ( lgamma(df / 2) - log(sqrt(df / 2)) - lgamma((df - 1) / 2) ))

# Approximate effect sizes using statistics (t)
print(t_to_d(t = c(-1.8188, 0.0981, 0.4514, -2.2114, -1.7010, 1.9169),
                   paired = F,
                   df_error = 114),
      digits = 3)

# Other ----
## GAMLj ----
# Effect sizes using the GAMLj
#devtools::install_github("gamlj/gamlj")
library(GAMLj3)

hexokinase_gam <- gamlj_lm(hexokinase ~ condition * timepoint + (1 | ID), 
                        data = metabolic_tidy,
                        emmeans = ~condition:timepoint,
                        posthoc = ~condition:timepoint,
                        adjust = "sidak",
                        posthoc_es = list('dm', 'g'),
                        d_ci = T)
hexokinase_gam$posthoc
hexokinase_gam$posthocEffectSize

## lme4 ----
# Using lme4
# Fit a linear mixed model with lme4 using condition and timepoint as both within factors
# and as their interaction.
library(lme4)
hexokinase_lmer <- lme4::lmer(hexokinase ~ condition * timepoint + (1 | ID), 
                              data = metabolic_tidy)

hexokinase_lmer2 <- 
  lme4::lmer(hexokinase ~ condition * timepoint + (1|ID/condition:timepoint), 
             data = hexokinase_dat,
             control=lmerControl(check.nobs.vs.nRE="ignore",
                                                                                                                        check.nobs.vs.nlev="ignore"),
             contrasts=list(condition=contr.sum, timepoint=contr.sum))

Anova(hexokinase_lmer2, type = "III", test = "F")

## car ----
# b. Using car package
library(car)

# Check the contrasts for each factor and their interaction.
contrasts(hexokinase_dat$condition)
contrasts(hexokinase_dat$timepoint)
contrasts(hexokinase_dat$condition:hexokinase_dat$timepoint)

# Fit a linear mixed model with car using condition and timepoint as both within factors
hexokinase_aovcar <- Anova(hexokinase_lm, 
                           idata = hexokinase_dat, 
                           idesign = ~ condition * timepoint,
                           multivariate = FALSE,
                           type = "III")

# This doesn't seem to provide the same results as the above two ways.

# agricolae package
library(agricolae)

# Use LSD.test() to get the pairwise comparisons for the interaction.
LSD.test(hexokinase_lm, 
         trt = c("condition", "timepoint"), 
         p.adj = "bonferroni",
         group = TRUE, 
         console = TRUE)

# Miscellaneous ----
## Base R ----
# Run many lm models
dv <- metabolic_tidy[, 4:6]
output <- lm(as.matrix(dv) ~ condition * timepoint, data = metabolic_tidy, 
             contrasts = list(condition = contr.sum, timepoint = contr.sum))
summary(aov(output))

rm(list = c('output','dv'))
rm(output, dv)

# Now fit anova's using aov()
dv <- metabolic_tidy[, 4:6]
output <- aov(as.matrix(dv) ~ condition * timepoint + Error(ID/(condition * timepoint)), 
              data = metabolic_tidy)
summary(aov(output))

# Now use the afex package
result <- lapply(names(metabolic_tidy)[4:6], function(x) { #ncol(metabolic_tidy)
  aov_mod <- aov_ez(id = "ID", 
                    dv = x, 
                    data = metabolic_tidy, 
                    within = c("condition", "timepoint"),
                    anova_table = list(es = "pes"))
   (emms1 <- emmeans(aov_mod, ~ condition * timepoint))
   update(pairs(emms1), adjust = "sidak")
})
result

metabolite_names <- colnames(metabolic_tidy)[4:ncol(metabolic_tidy)]

# Could not make the ez package to work.

# Now with the tidyverse
## Tidyverse ----

# Tidy the data
# Sanity check 17*20*2*3 = 2040
metabolic_tidier <- 
  pivot_longer(data = metabolic_tidy,
               cols = !c(ID, condition, timepoint),
               cols_vary = "fastest",
               values_to = "Measurement",
               names_to = "Molecule")

# Use the afex package
result <- metabolic_tidier %>%
  group_by(Molecule) %>%
  nest() %>%
  mutate(aov = map(data, ~ aov_ez(id = "ID", 
                                  dv = 'Measurement', 
                                  data = .x, 
                                  within = c("condition", "timepoint"),
                                  anova_table = list(es = "pes"))))
result$aov

# Use the ezANOVA package
result <- metabolic_tidier %>%
  group_by(Molecule) %>%
  nest() %>%
  mutate(aov = map(data, ~ ezANOVA(data = .x, 
                                   dv = Measurement, 
                                   wid = ID, 
                                   within = c(condition, timepoint), 
                                   detailed = TRUE, 
                                   type = 3)))
result$aov

# Change the name of the lists
metabolite_names <- colnames(metabolic_tidy)[4:ncol(metabolic_tidy)]
names(result$aov) <- metabolite_names
result$aov
