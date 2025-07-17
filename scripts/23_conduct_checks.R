
source("scripts/00_config.R")

library(tidyverse)
library(bayesplot)
library(tidybayes)
library(posterior)
library(rstan)

# Import data (variable names, model outputs)
a_imp_df <- colnames(readRDS("mars/conduct_model/latest_stan_conduct.rds")$X)
a_imp_df <- a_imp_df[1:19]
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

# Plot betas for rate
betas <- a_imp_res %>% 
  filter(grepl("^beta\\[", variable)) %>% 
  separate(variable, into = c("beta", "order"), sep = "\\[") %>% 
  mutate(order = gsub("\\]", "", order)) %>% 
  cbind(data.frame(var = a_imp_df)) %>% 
  select(var, everything())

betas_plot <- betas %>% 
  ggplot(aes(x = mean, xmin = `2.5%`, xmax = `97.5%`, y = fct_rev(var))) +
  geom_point() +
  geom_linerange() +
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
  scale_x_continuous(n.breaks = 15) +
  theme_bw() +
  labs(x = "Mean point estimate (95% credible interval", y = NULL)

betas_plot

ggsave(
  "mars/conduct_model/outputs/rate_betas.png",
  betas_plot,
  width = 8,
  height = 8,
  units = "in"
)
