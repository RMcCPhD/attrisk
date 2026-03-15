
# Summary statistics and plots
# Conduct-stage model

source("scripts/00_config.R")
library(tidyverse)
library(posterior)
library(bayesplot)
library(flexsurv)
library(metafor)
library(broom)

# Import trial-level Cox model outputs
a_imp_res <- read_csv("vivli/01_res.csv")
a_imp_diag <- read_csv("vivli/02_diag.csv")
a_imp_vcov <- read_csv("vivli/03_cor.csv")

# Tidy covar names
z_covars <- c(
  "Age (in decades)", "BMI (kg/m2)", "HDL cholesterol (mmol/L)",
  "Total cholesterol (mmol/L)", "N comorbidities",
  "Duration of T2DM (years)", "eGFR (ml/min/1.73m2)", "HbA1c (%)",
  "N conc. medication classes", "Asian", "Black/African American",
  "Other race/ethnicity", "Systolic BP (mm/Hg)", "Male", "Current smoker",
  "Treatment: Active comparator"
)

z_covars_lvl <- c(
  "Treatment: Active comparator",
  "Age (in decades)", "Male", "Asian", "Black/African American",
  "Other race/ethnicity", "Current smoker", "BMI (kg/m2)",
  "Systolic BP (mm/Hg)", "HDL cholesterol (mmol/L)",
  "Total cholesterol (mmol/L)",
  "Duration of T2DM (years)", "HbA1c (%)", "eGFR (ml/min/1.73m2)",
  "N conc. medication classes", "N comorbidities"
)

b_covars <- a_imp_res %>% 
  count(term) %>% 
  cbind(z_covars) %>% 
  as_tibble() %>% 
  select(-n) %>% 
  mutate(
    z_covars = factor(
      z_covars,
      levels = z_covars_lvl
    )
  ) %>%
  rename(covars = z_covars) %>%
  arrange(covars)

# IPD model: Trial-level log-hazard ratios -------------------------------------

# Convergences issues with covariates or missing category levels results in
# an incomplete set of covariates across trials
# Addressed by inserting via expansion grid, setting estimates to 0 and
# blown up SEs
# In MA, VCOVs blown up to promote shrinkage towards zero

# Tidy point estimates
b_res_tidy <- b_covars %>% 
  left_join(a_imp_res) %>%
  select(nct_id, covars, estimate, std.error, statistic, p.value) %>% 
  arrange(nct_id, covars)

# Get unique model terms
c_unq <- unique(b_res_tidy$covars)

# Complete terms set (i.e. were missing from IPD, convg issues)
c_covars_grid <- expand.grid(
  nct_id = unique(b_res_tidy$nct_id),
  covars = c_unq
  ) %>% 
  as_tibble()

# Join grid to estimates
# Add zero estimates and large SEs where missing
c_compl <- b_res_tidy %>%
  full_join(c_covars_grid) %>% 
  mutate(
    estimate = if_else(is.na(estimate), 0, estimate),
    std.error = if_else(is.na(std.error), 99999, std.error)
  ) %>% 
  arrange(nct_id, covars)

# Round to 3DP
# Calculate confidence intervals (set to NA where covariate absent from model)
c_rdy <- c_compl %>% 
  mutate(
    lci = estimate - 1.96 * std.error,
    uci = estimate + 1.96 * std.error,
    across(lci:uci, ~ if_else(std.error == 99999, NA, std.error)),
    across(estimate:uci, ~ round(., 3)),
  ) %>% 
  select(
    nct_id, covars, estimate, std.error, lci, uci,
    everything()
  )

# Save
write_csv(c_rdy, "outputs/summaries/conduct_model/ipd_cox_model_res.csv")

# Summarise model diagnostics
c_diag <- a_imp_diag
write_csv(c_diag, "outputs/summaries/conduct_model/ipd_cox_model_diag.csv")

# MA: Pooled log-hazard and hazard ratios --------------------------------------

# Correctly order estimates
d_lvl <- a_imp_res %>% count(term) %>% pull(term)

d_est <- c_rdy %>% 
  left_join(b_covars) %>% 
  mutate(term = factor(term, levels = c(d_lvl))) %>% 
  arrange(nct_id, term)

# Get VCOV
d_vcov <- readRDS("processed_data/conduct_cox_vcov.rds")

# Fit meta-regression with random effects for covariates and trials
# Random effects for term and trial
d_ma_id <- d_est$nct_id
d_ma_params <- d_est$term
d_ma_ests <- d_est$estimate
d_ma_vcov <- as.matrix(bdiag(d_vcov$cov_matrix))

d_ma <- rma.mv(
  yi = d_ma_ests,
  V = d_ma_vcov,
  mods = ~ factor(d_ma_params) - 1,
  random = ~ d_ma_params | d_ma_id
)

# Checks
summary(d_ma)
d_ma$tau2
d_ma$rho

# Tidy results
d_ma_tidy <- tidy(d_ma) %>% 
  mutate(term = gsub("factor\\(d_ma_params\\)", "", term)) %>% 
  mutate(
    lo = estimate - 1.96 * std.error,
    hi = estimate + 1.96 * std.error,
    hr = exp(estimate),
    lo_hr = exp(lo),
    hi_hr = exp(hi)
  ) %>% 
  select(-type) %>% 
  left_join(b_covars) %>% 
  arrange(covars) %>% 
  select(
    covars, estimate, lo, hi, hr, lo_hr, hi_hr,
    std.error, statistic, p.value
  ) %>% 
  rename(
    loghr = estimate,
    loghr_lo = lo,
    loghr_hi = hi
  ) %>% 
  mutate(
    across(loghr:statistic, ~ round(., 2)),
    p.value = round(p.value, 3)
  )

# Plot
d_ma_plot <- d_ma_tidy %>% 
  ggplot(aes(x = hr, xmin = lo_hr, xmax = hi_hr, y = fct_rev(covars))) +
  geom_point() +
  geom_linerange() +
  geom_vline(xintercept = 1, colour = "red", linetype = "dashed") +
  scale_x_continuous(limits = c(0.4, 1.6), n.breaks = 7) +
  theme_classic() +
  labs(
    x = "Mean effect on attrition risk (95% confidence intervals)", 
    y = NULL
  )

# Save
write_csv(d_ma_tidy, "outputs/summaries/conduct_model/ipd_ma_estimates.csv")

ggsave(
  "outputs/plots/conduct_ipd_chars_risk.png",
  d_ma_plot,
  height = 6,
  width = 6,
  units = "in"
)

# Calibration-in-the-large -----------------------------------------------------

# Import materials from CITL
e_citl_draws <- readRDS("outputs/conduct_model/citl/delta_draws.rds")

# Import CITL summary
# Add trial ID
e_citl_sum <- read_csv("outputs/conduct_model/citl/delta_risk_sum.csv") %>% 
  cbind(
    tibble(
      trial_id = readRDS("processed_data/latest_dataset.rds") %>% 
        pull(trial_id)
    )
  ) %>% 
  as_tibble()

# Summarise distribution of delta (intercept shifts)
e_citl_overall <- tibble(
  mean   = mean(e_citl_draws),
  median = median(e_citl_draws),
  sd     = sd(e_citl_draws),
  lo     = quantile(e_citl_draws, 0.025),
  hi     = quantile(e_citl_draws, 0.975)
)

# Summarise trial-level delta draws
e_citl_trial <- tibble(
  trial_id = e_citl_sum$trial_id,
  trial_arm = seq_len(ncol(e_citl_draws)),
  delta_mean = colMeans(e_citl_draws),
  delta_median = apply(e_citl_draws, 2, median),
  delta_lo = apply(e_citl_draws, 2, quantile, 0.025),
  delta_hi = apply(e_citl_draws, 2, quantile, 0.975)
)

# Summarise distribution of calibration shifts
e_citl_dist <- tibble(delta = as.vector(e_citl_draws)) %>%
  ggplot(aes(x = delta)) +
  geom_density(fill = "steelblue", alpha = .4) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic() +
  labs(
    x = "Calibration intercept shift (δ)",
    y = "Density"
  )

# Summarise pre and post-CITL predicted risk
e_citl_risk <- e_citl_sum %>% 
  reframe(
    obs_mean  = mean(p_obs),
    pred_pre_mean  = mean(p_pred),
    pred_post_mean = mean(p_citl)
  )

# Summarise absolute calibration error
e_citl_error <- e_citl_sum %>%
  mutate(
    err_pre  = abs(p_pred - p_obs),
    err_post = abs(p_citl - p_obs)
  ) %>%
  reframe(
    mae_pre  = mean(err_pre),
    mae_post = mean(err_post)
  )

# Calibration plot
e_citl_plot_dat <- tibble(
  arm = rep(e_citl_sum$arm, 2),
  obs = rep(e_citl_sum$p_obs, 2),
  pred = c(e_citl_sum$p_pred, e_citl_sum$p_citl),
  stage = rep(c("Pre-CITL", "Post-CITL"), each = nrow(e_citl_sum))
)

e_citl_plot <- e_citl_plot_dat %>%
  ggplot(aes(x = obs, y = pred, colour = stage)) +
  geom_abline(linetype = 2) +
  geom_point(size = 0.5) +
  theme_classic() +
  labs(
    x = "Observed attrition risk",
    y = "Predicted attrition risk",
    colour = NULL
  ) +
  scale_x_continuous(limits = c(0, 1), n.breaks = 5)

# Save summaries
write_csv(e_citl_overall, "outputs/summaries/conduct_model/citl_delta_sum_overall.csv")
write_csv(e_citl_sum, "outputs/summaries/conduct_model/citl_delta_sum_trial.csv")
write_csv(e_citl_risk, "outputs/summaries/conduct_model/citl_risk_sum.csv")
write_csv(e_citl_error, "outputs/summaries/conduct_model/citl_error.csv")

write_csv(e_citl_trial, "outputs/summaries/conduct_model/citl_tmp.csv")

# Save plots
ggsave(
  "outputs/plots/conduct_citl_delta.png",
  e_citl_dist,
  width = 4,
  height = 4,
  units = "in"
)

ggsave(
  "outputs/plots/conduct_citl_calibration.png",
  e_citl_plot,
  width = 4,
  height = 4,
  units = "in"
)

# Case scenario attrition risk prediction --------------------------------------
