
rm(list = ls())
gc()

# Packages
library(tidyverse)
library(rstan)
library(bayesplot)

# Config
options(scipen = 999)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Import stan data
a_imp_dat <- readRDS("stan_dat.rds")

# Compile stan model
b_compile <- stan_model("mdl.stan")

# Fit model
c_fit <- sampling(
  b_compile,
  data = a_imp_dat,
  warmup = 1000,
  iter = 3000,
  chains = 4,
  seed = 123,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

saveRDS(c_fit, "design_fit.rds")
