
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

# Import stan data
a_imp_dat <- readRDS("stan_dat.rds")

# Compile stan model
b_compile <- stan_model("mdl_int.stan")

# Fit model
c_fit <- sampling(
  b_compile,
  data = a_imp_dat,
  warmup = 1000,
  iter = 4000,
  chains = 4,
  seed = 123,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

saveRDS(c_fit, "int_fit.rds")

# Summarise
d_draws <- as_draws_df(c_fit)
d_sum_draws <- summarise_draws(d_draws)

# Check for divergences
check_hmc_diagnostics(c_fit) # Nothing abnormal, i.e. 0 divergences

# Split draws into a separate column
d_sum <- d_sum_draws %>% 
  mutate(variable = gsub("\\[[0-9,]+\\]", "", variable)) %>% 
  group_by(variable) %>% 
  mutate(draw = case_when(n() == 1 ~ NA, TRUE ~ row_number())) %>% 
  select(variable, draw, everything()) %>% 
  ungroup()

d_sum_conv <- d_sum %>% 
  group_by(variable) %>% 
  reframe(
    n = n(),
    mean_rhat = mean(rhat, na.rm = TRUE),
    mean_ess = mean(ess_bulk, na.rm = TRUE)
  )

write_csv(d_sum_conv, "sum_rhat_ess.csv")

# Extract posterior estimates
e_sum <- summarise_draws(
  d_draws,
  mean, median, sd,
  ~ quantile(.x, probs = c(0.025, 0.975))
)

saveRDS(e_sum, "int_sum_draws.rds")

# Examine r recovery
ppc_dens <- ppc_dens_overlay(a_imp_dat$r, as.matrix(c_fit, pars = "r_new")[1:361, ])
ggsave("r_recovery_overall.jpg", ppc_dens, width = 8, height = 8, units = "in")

ppc_int <- ppc_intervals(a_imp_dat$r, as.matrix(c_fit, pars = "r_new")[1:361, ])
ggsave("r_recovery_ind.jpg", ppc_int, width = 8, height = 8, units = "in")
