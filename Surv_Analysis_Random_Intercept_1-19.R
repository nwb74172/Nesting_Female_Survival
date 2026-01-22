###############################################################################L


# Female turkey daily survival rate analysis; Random effect of site


# Code written by: Allison Keever
#                 akeever@tntech.edu

# Data cleaned and supplied by: Nicholas Bakner
#                               nbakner@tntech.edu

# Last updated: 01/19/2026


###############################################################################L


# Prep work space---------------------------------------------------------------

# Load packages
library(tidyverse)
library(jagsUI)
library(tidybayes)

# Load data
datum <- read.csv("Data/Veg_Survival_Encounter_History_Site.csv") %>%
  select(-X)


# Source model code
source("Rscripts/DSR_Model.R")


# Pre-analysis: data summaries and correlations---------------------------------

# Look at ages in data. Because there are so few juveniles, adults will be in 
# the intercept. Also, there are 4 unknown birds, and we will just try to 
# impute the ages of those. If that fails, we will toss them. 
table(datum$age)


# Look at correlation chart. First have to turn age into numeric (A = 0, J = 1).
# Average veg and VO veg were correlated (r = 0.8), so average veg was dropped.
# This figure was saved. 
PerformanceAnalytics::chart.Correlation(
  datum %>%
    mutate(age = case_when(age == "A" ~ 0,
                           age == "J" ~ 1, 
                           TRUE ~ NA)) %>%
    select(gc, cc, tph, maxveg, avgveg, voveg, age, year))


# Create summary table of covariates
datum %>% 
  summarise(mean_gc = mean(gc, na.rm = TRUE), 
            sd_gc = sd(gc, na.rm = TRUE),
            mean_cc = mean(cc, na.rm = TRUE), 
            sd_cc = sd(cc, na.rm = TRUE),
            mean_tph = mean(tph, na.rm = TRUE), 
            sd_tph = sd(tph, na.rm = TRUE),
            mean_maxveg = mean(maxveg, na.rm = TRUE), 
            sd_maxveg = sd(maxveg, na.rm = TRUE),
            mean_voveg = mean(voveg, na.rm = TRUE), 
            sd_voveg = sd(voveg, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), 
               names_to = c("stat", "Variable"), 
               names_sep = "_",
               values_to = "Value") %>% 
  pivot_wider(names_from = "stat", 
              values_from = "Value")


# Run analysis------------------------------------------------------------------

# Find number of observations per female (days alive and nesting)
nobs <- apply(datum %>% 
                select(starts_with("V"), -voveg), 
              1, function(x) sum(!is.na(x)))

# Set up jags data
jags_data <- list(
  "ninds" = nrow(datum), "nobs" = nobs, "nsites" = length(unique(datum$site)),
  "y" = datum %>% 
    select(starts_with("V"), -voveg),
  "Age" = datum %>% 
    transmute(age = case_when(age == "A" ~ 0,
                              age == "J" ~ 1, 
                              TRUE ~ NA_integer_)), 
  "Site" = datum %>% 
    mutate(Site = case_when(site == "BFG" ~ 1, 
                            site == "CC" ~ 2, 
                            site == "SELA" ~ 3, 
                            site == "SL" ~ 4, 
                            site == "SRS" ~ 5, 
                            site == "WEBB" ~ 6, 
                            site == "WLA" ~ 7, 
                            TRUE ~ 8)) %>% 
    pull(Site),
  "CC" = datum$cc, "TPH" = datum$tph, "Max" = datum$maxveg, "VO" = datum$voveg, 
  "thr" = seq(0, 1, 0.005), "GC" = datum$gc, 
  "nzeros" = length(which(datum %>% 
                            select(starts_with("V"), -voveg)==0)))


# Initial values
inits <- function(){list()}


# Parameters to track
params <- c("b0", "b.juv", "b.cc", "b.gc", "b.tph", "b.maxveg", "b.voveg", 
            "cuml.q", "p", "auc", "sd.site", "loglik")


# MCMC settings
ni <- 5000
nmax <- 200000
nb <- 20000
nc <- 3
nt <- 5


# Run model 
out <- autojags(data = jags_data, inits = inits, parameters.to.save = params, 
            model.file = "DSR_Female_Random_Intercept.txt", n.chains = nc, 
            iter.increment = ni, n.burnin = nb, n.thin = nt, 
            max.iter = nmax, parallel = TRUE)

# Look at betas and CRIs 

df <- out %>% 
  gather_draws(b0, b.juv, b.cc, b.gc, b.tph, b.maxveg, b.voveg,  
               sd.site) %>%
  mean_qi(.width = 0.95)

# Determine if ground cover is biologically relevant

# Create predictive frame

# Covariate means (use na.rm=TRUE)
mu_cc     <- mean(datum$cc, na.rm = TRUE)
mu_tph    <- mean(datum$tph, na.rm = TRUE)
mu_maxveg <- mean(datum$maxveg, na.rm = TRUE)
mu_voveg  <- mean(datum$voveg, na.rm = TRUE)

# Set age scenario (adult=0; juvenile=1)
age_scn <- 1  # adult

# GC sequence across the observed (central) range
gc_seq <- seq(
  quantile(datum$gc, 0.05, na.rm = TRUE),
  quantile(datum$gc, 0.95, na.rm = TRUE),
  length.out = 150
)

pred_grid <- expand_grid(
  gc = gc_seq,
  cc = mu_cc,
  tph = mu_tph,
  maxveg = mu_maxveg,
  voveg = mu_voveg,
  age = age_scn
)

# Extract posterior draws as a data frame
draws <- as_draws_df(out$samples)

# Keep only the fixed effects you need
beta <- draws %>%
  transmute(
    .draw = row_number(),
    b0 = b0,
    b_gc = b.gc,
    b_cc = b.cc,
    b_tph = b.tph,
    b_maxveg = b.maxveg,
    b_voveg = b.voveg,
    b_juv = b.juv
  )

# Predict daily survival and 33-day survival across GC grid
pred_draws <- pred_grid %>%
  crossing(beta) %>%
  mutate(
    eta = b0 +
      b_juv   * age +
      b_cc    * cc +
      b_tph   * tph +
      b_maxveg* maxveg +
      b_voveg * voveg +
      b_gc    * gc,
    q_daily = plogis(eta),
    s33 = q_daily^33
  )

pred_summ <- pred_draws %>%
  group_by(gc) %>%
  median_qi(q_daily, s33, .width = 0.90)

# Plot 33-day survival vs GC
ggplot(pred_summ, aes(x = gc, y = s33)) +
  geom_ribbon(aes(ymin = s33.lower, ymax = s33.upper), alpha = 0.25) +
  geom_line(linewidth = 1) +
  labs(
    x = "Ground cover (%)",
    y = "Predicted incubation survival",
  ) +
  theme_classic()
