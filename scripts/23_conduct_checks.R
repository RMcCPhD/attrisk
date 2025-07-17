

source("scripts/00_config.R")

library(tidyverse)
library(bayesplot)
library(tidybayes)
library(posterior)
library(rstan)

a_imp_df <- colnames(readRDS("mars/conduct_model/latest_stan_conduct.rds")$X)
a_imp_res <- readRDS("mars/conduct_model/outputs/conduct_sum_draws.rds")
a_imp_mdl <- readRDS("mars/conduct_model/outputs/conduct_fit.rds")
a_imp_rhat <- read_csv("mars/conduct_model/outputs/sum_rhat_ess.csv")

# Trace plots for global parameters
trace_global <- traceplot(
  a_imp_mdl, 
  pars = c(
    "lp__", "mean_rate", "sd_rate", "mu_mu1", "tau_mu1", "tau_beta"
  )
)

ggsave(
  "mars/conduct_model/outputs/trace_global.png", 
  trace_global, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots for covariate betas
trace_beta <- traceplot(
  a_imp_mdl, 
  pars = "beta"
)

ggsave(
  "mars/conduct_model/outputs/trace_beta.png", 
  trace_beta, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots - first 56 of shape and rate
trace_a <- traceplot(a_imp_mdl, pars = paste0("a[", 1:56, "]"))
trace_b <- traceplot(a_imp_mdl, pars = paste0("b[", 1:56, "]"))

ggsave(
  "mars/conduct_model/outputs/trace_a.png", 
  trace_a, 
  width = 15, 
  height = 10,
  units = "in"
)

ggsave(
  "mars/conduct_model/outputs/trace_b.png", 
  trace_b, 
  width = 15, 
  height = 10,
  units = "in"
)
