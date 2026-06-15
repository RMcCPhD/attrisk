
# Test out the Stan model to identify any glaring
# model specification or data-related issues

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

# Import AgD
imp_dat <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Import Stan data
imp_mdl_dat <- readRDS("processed_data/stan_data.rds")

# Compile model
mdl_compile <- stan_model("scripts/04a_design_model.stan")

# Take first 30 trials
dat_tst <- imp_mdl_dat
dat_tst$I <- 30
dat_tst$r <- dat_tst$r[seq_len(30)]
dat_tst$N <- dat_tst$N[seq_len(30)]
dat_tst$dt <- dat_tst$dt[seq_len(30)]
dat_tst$X <- dat_tst$X[seq_len(30), , drop = FALSE]

# Fit model
# Converged, no warnings about divergence or ESS
mdl_fit <- sampling(
  mdl_compile,
  data = dat_tst,
  warmup = 1000,
  iter = 4000,
  chains = 4,
  seed = 123,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

# Get posterior summaries
ppc_out <- extract(
  mdl_fit, 
  pars = c("r_new", "p", "log_lik", "beta")
)

ppc_tbl <- tibble(
  row = seq_len(dat_tst$I),
  r_obs = dat_tst$r,
  N = dat_tst$N,
  obs_risk = r_obs / N,
  r_pred_mean = colMeans(ppc_out$r_new),
  r_pred_q025 = apply(ppc_out$r_new, 2, quantile, 0.025),
  r_pred_q50 = apply(ppc_out$r_new, 2, quantile, 0.50),
  r_pred_q975 = apply(ppc_out$r_new, 2, quantile, 0.975),
  p_mean = colMeans(ppc_out$p),
  p_q025 = apply(ppc_out$p, 2, quantile, 0.025),
  p_q50 = apply(ppc_out$p, 2, quantile, 0.50),
  p_q975 = apply(ppc_out$p, 2, quantile, 0.975)
)

write_csv(
  ppc_tbl %>% mutate(across(everything(), ~ round(., 3))), 
  "outputs/design/model_test/posterior_summary.csv"
)

# Performance summary
# Model recovers observed counts very well
ppc_performance <- ppc_tbl %>%
  mutate(
    covered_95 = r_obs >= r_pred_q025 & r_obs <= r_pred_q975,
    error = r_pred_mean - r_obs,
    abs_error = abs(error),
    sq_error = error^2,
    rel_error = (r_pred_mean - r_obs) / pmax(r_obs, 1),
    abs_rel_error = abs(rel_error)
  ) %>%
  reframe(
    n = n(),
    covered_95 = mean(covered_95),
    mean_error = mean(error),
    rel_error = mean(rel_error),
    mae = mean(abs_error),
    rmse = sqrt(mean(sq_error)),
    cor_obs_pred = cor(r_obs, r_pred_mean)
  ) %>% 
  mutate(across(everything(), ~ round(., 3)))

write_csv(ppc_performance, "outputs/design/model_test/posterior_performance.csv")

# Plot posterior recovery
# Confirms adequate recovery
ppc_plot <- ppc_tbl %>% 
  ggplot(aes(x = row, y = r_obs)) +
  geom_point(alpha = 0.7) +
  geom_point(
    aes(x = row, y = r_pred_mean), 
    colour = "skyblue",
  ) +
  geom_errorbar(
    aes(ymin = r_pred_q025, ymax = r_pred_q975), 
    colour = "skyblue",
    width = 0.1
  ) +
  theme_classic() +
  labs(
    x = "Trial arm", 
    y = "Attrition count",
    title = "Observed versus posterior predicted attrition counts"
  )

ggsave(
  "outputs/design/model_test/posterior_count_recovery.png",
  ppc_plot,
  width = 10,
  height = 4,
  units = "in"
)

# Get posterior covariate effect estimates
res_names <- colnames(readRDS("processed_data/agd_covar_mat.rds"))
res_covars <- as.data.frame(ppc_out$beta)
names(res_covars) <- res_names

# Summarise covariate effect estimates
# Subset of trials, don't expect findings to be telling
res_sum <- res_covars %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  reframe(
    mu = mean(value),
    q025 = quantile(value, 0.025),
    q975 = quantile(value, 0.975)
  )

res_plot <- res_sum %>% 
  ggplot(aes(x = mu, y = fct_rev(name))) +
  geom_point() +
  geom_errorbar(aes(xmin = q025, xmax = q975)) +
  theme_classic() +
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
  labs(
    x = "Trial arm", 
    y = "Effect on log-rate",
  )

ggsave(
  "outputs/design/model_test/posterior_effects_log_rate.png",
  res_plot,
  width = 8,
  height = 6,
  units = "in"
)
