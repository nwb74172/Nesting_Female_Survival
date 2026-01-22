###############################################################################L


# Female turkey daily survival rate analysis


# Code written by: Allison Keever
#                 akeever@tntech.edu

# Data cleaned and supplied by: Nick Bakner
#                               nbakner@tntech.edu

# Last updated: 01/08/2026


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
            "cuml.q", "p", "auc", "u.int", "u.gc", "u.cc", "u.tph", "u.maxveg",
            "u.voveg", "sd.int", "sd.gc", "sd.cc", "sd.tph", "sd.maxveg",
            "sd.voveg", "loglik")


# MCMC settings
ni <- 5000
nmax <- 200000
nb <- 20000
nc <- 3
nt <- 5


# Run model 
out <- autojags(data = jags_data, inits = inits, parameters.to.save = params, 
            model.file = "DSR_Female_Random_Slope.txt", n.chains = nc, 
            iter.increment = ni, n.burnin = nb, n.thin = nt, 
            max.iter = nmax, parallel = TRUE)

# Look at betas and CRIs 

df <- out %>% 
  gather_draws(b0, b.juv, b.cc, b.gc, b.tph, b.maxveg, b.voveg,  
               u.int[1:7], u.gc[1:7], u.cc[1:7], u.tph[1:7], 
               u.maxveg[1:7], u.voveg[1:7], sd.int, sd.gc,
               sd.cc, sd.tph, sd.maxveg, sd.voveg) %>%
  mean_qi(.width = 0.95)

df <- as.data.frame(out$summary)
