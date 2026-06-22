
# Posterior predictive checks for full design stage model
# Trace plots not generated because Rhat = 1 achieved for all parameters
# Assessing:
# - Recovery of attrition counts
# - Predictive performance
# - Covariate effect estimates (not PPC but here for initial check)

rm(list = ls())
gc()

# Packages
library(tidyverse)
library(rstan)
library(bayesplot)
library(posterior)

# Import Stan data
imp_dat <- readRDS("processed_data/stan_data.rds")

# Import fitted model
imp_mdl <- readRDS("outputs/design/design_model.rds")

# Get posterior summaries
ppc_out <- extract(
  imp_mdl, 
  pars = c("r_new", "p", "log_lik", "beta")
)

ppc_tbl <- tibble(
  row = seq_len(imp_dat$I),
  r_obs = imp_dat$r,
  N = imp_dat$N,
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
  "outputs/design/posterior_summary.csv"
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

write_csv(ppc_performance, "outputs/design/posterior_performance.csv")

# Plot posterior recovery
# Confirms adequate recovery
ppc_plot <- ppc_tbl %>% 
  ggplot(aes(x = row, y = r_obs)) +
  geom_point(size = 0.5) +
  geom_point(
    aes(x = row, y = r_pred_mean), 
    colour = "skyblue",
    size = 0.5,
    alpha = 0.7
  ) +
  geom_linerange(
    aes(ymin = r_pred_q025, ymax = r_pred_q975), 
    colour = "skyblue",
    alpha = 0.7
  ) +
  theme_classic() +
  labs(
    x = "Trial arm", 
    y = "Attrition count",
    # title = "Observed versus posterior predicted attrition counts"
  )

ggsave(
  "outputs/design/posterior_count_recovery.png",
  ppc_plot,
  width = 10,
  height = 4,
  units = "in"
)

# Get tidy covariate names
res_names <- colnames(readRDS("processed_data/agd_covar_mat.rds"))
# write_csv(
#   tibble(raw_name = res_names), 
#   "scratch_data/data_processing/tidy_covariate_names.csv"
# )

# Get posterior covariate effect estimates
res_covars <- as.data.frame(ppc_out$beta)
names(res_covars) <- res_names
z_names <- "scratch_data/data_processing/tidy_covariate_names.csv"

res_names <- read_csv(z_names) %>% 
  mutate(
    tidy_name = factor(tidy_name, levels = unique(tidy_name))
  )

# Summarise covariate effect estimates
# Subset of trials, don't expect findings to be telling
res_sum <- res_covars %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  reframe(
    mu = mean(value),
    q025 = quantile(value, 0.025),
    q975 = quantile(value, 0.975)
  ) %>% 
  left_join(res_names %>% rename(name = raw_name)) %>% 
  arrange(tidy_name)

res_plot <- res_sum %>% 
  ggplot(aes(x = mu, y = fct_rev(tidy_name))) +
  geom_point() +
  geom_linerange(aes(xmin = q025, xmax = q975)) +
  scale_x_continuous(limits = c(-1.2, 1.2), n.breaks = 6) +
  theme_classic() +
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
  labs(
    x = "Mean effect on Gompertz log-rate (95% credible intervals)", 
    y = NULL,
  )

ggsave(
  "outputs/design/posterior_effects_log_rate.png",
  res_plot,
  width = 10,
  height = 6,
  units = "in"
)
