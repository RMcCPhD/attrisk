
# Fit final design-stage model
# Script is run on a high performance computing unit as runtime on local
# machine would be significantly longer
# Can run locally but takes an extremely long time due to computing demands

rm(list = ls())
gc()

# Packages
library(tidyverse)
library(rstan)
library(bayesplot)
library(posterior)

# Config
options(scipen = 999)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Import Stan data
imp_mdl_dat <- readRDS("stan_data.rds")

# Compile model
mdl_compile <- stan_model("04a_design_model.stan")

# Fit model
mdl_fit <- sampling(
  mdl_compile,
  data = imp_mdl_dat,
  warmup = 1000,
  iter = 4000,
  chains = 4,
  seed = 123,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

saveRDS(mdl_fit, "design_model.rds")