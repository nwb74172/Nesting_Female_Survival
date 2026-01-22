###############################################################################L


# Predictive Performance of random effect of site and random slopes of site

# Data cleaned and supplied by: Nicholas Bakner
#                               nbakner@tntech.edu

# Last updated: 01/19/2026


###############################################################################L

# Load packages
library(tidyverse)
library(jagsUI)
library(tidybayes)
library(loo)
library(posterior)

# Data used to create matrices for comparison are from out in 
# Surv_Analysis_Random_Intercept and Surv_Analysis_Random_Slope

# Preform WAIC to compare models
extract_loglik_matrix <- function(jags_out) {
  
  # Convert to draws_df
  draws <- as_draws_df(jags_out$samples)
  
  # Pull loglik columns only
  loglik_draws <- draws %>%
    select(starts_with("loglik["))
  
  # Convert to matrix: iterations x observations
  loglik_mat <- as.matrix(loglik_draws)
  
  return(loglik_mat)
}

# Calculate WAIC intercept
loglik_int <- extract_loglik_matrix(out)

waic_int <- waic(loglik_int)
print(waic_int)

# Remove to put in the slope data (accidently named things the same)
rm(list = setdiff(ls(), c("loglik_int", "loglik_slope", "waic_slope", 
                          "waic_int")))

# Calculate WAIC for slope
loglik_slope <- extract_loglik_matrix(out)

waic_slope <- waic(loglik_slope)
print(waic_slope)

# Comparison
loo_compare(waic_int, waic_slope)

################################################################################
# Conclusion: Including random slopes does not increase the predictive 
# performance of the model. Therefore, we go with the least complex model
# which includes site as a random effect.
################################################################################