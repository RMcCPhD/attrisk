
# Summary statistics and plots
# Design stage model

source("scripts/00_config.R")
library(tidyverse)
library(posterior)
library(bayesplot)
library(flexsurv)

# Import AgD, posterior draws and model
a_imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]
a_imp_fit <- readRDS("outputs/design/design_model.rds")

# Import STAN dataset to identify covariates
a_imp_stan <- readRDS("processed_data/stan_data.rds")

# Get covariate names
b_var_names <- tibble(
  covar = colnames(a_imp_stan$X),
  draw = seq(colnames(a_imp_stan$X))
)

# Add tidy names manually and import
# write_csv(b_var_names, "scratch_data/tidy_covar_names.csv")

b_var_names <- read_csv("scratch_data/data_processing/tidy_covariate_names.csv") %>% 
  filter(!raw_name == "upper_age80") %>% 
  mutate(
    tidy_name = factor(tidy_name, levels = unique(tidy_name)),
    variable = paste0("beta[", row_number(), "]")
  )

# General posterior summaries --------------------------------------------------

# Summarise draws
b_draws_array <- as_draws_df(a_imp_fit)
b_sum_draws <- summarise_draws(
  b_draws_array,
  mean, median, sd, rhat, ess_bulk, ess_tail,
  ~ quantile(.x, probs = c(0.025, 0.975))
)

# Tidy draws summary
b_sum_tidy <- b_sum_draws %>% 
  mutate(parameter = gsub("\\[[0-9,]+\\]", "", variable)) %>% 
  group_by(parameter) %>% 
  mutate(draw = case_when(n() == 1 ~ NA, TRUE ~ row_number())) %>% 
  select(parameter, draw, everything()) %>% 
  ungroup()

# Add covariate names
b_add_names <- b_sum_tidy %>% 
  filter(!parameter %in% c("beta", "lambda_beta", "z_beta")) %>% 
  full_join(
    b_sum_tidy %>% 
      filter(parameter %in% c("beta", "lambda_beta", "z_beta")) %>% 
      left_join(b_var_names)
  ) %>% 
  arrange(parameter, draw) %>% 
  select(parameter, draw, tidy_name, everything())

# Keep core parameters
b_sum_core <- b_add_names %>% 
  select(-variable) %>% 
  filter(
    parameter %in% c("a", "b", "log_b", "beta", "p", "r_new")
  )

# Posterior summaries: Trial characteristics -----------------------------------

# Isolate characteristics
c_char_sum <- b_sum_core %>% 
  filter(parameter == "beta") %>% 
  select(-c(parameter, draw)) %>% 
  mutate(across(mean:`97.5%`, ~ round(., 2)))

# Plot
c_char_plot <- c_char_sum %>% 
  ggplot(aes(x = mean, xmin = `2.5%`, xmax = `97.5%`, y = fct_rev(covar_name))) +
  geom_point() +
  geom_linerange() +
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
  scale_x_continuous(limits = c(-0.8, 0.8), n.breaks = 9) +
  theme_classic() +
  labs(
    x = "Mean effect on Gompertz log-rate (95% credible intervals)", 
    y = NULL
  )

# Save
write_csv(c_char_sum, "outputs/manuscript_summaries/design_trial_chars.csv")

# Posterior summaries: Model parameters ----------------------------------------

# Isolate model parameters
d_mdl_sum <- b_add_names %>% 
  select(-variable) %>% 
  filter(!parameter == "beta") %>% 
  select(-tidy_name, -raw_name)

# Add nct_id (1:1 order with STAN data)
# Split by type
d_mdl_nest <- d_mdl_sum %>% 
  left_join(
    a_imp_agd %>% 
      select(trial_id) %>% 
      distinct() %>% 
      mutate(draw = seq(1:length(trial_id)))
  ) %>% 
  select(-draw) %>% 
  group_by(trial_id, parameter) %>% 
  nest(.key = "draws")

# Check parameters
d_mdl_nest %>% ungroup() %>% count(parameter)

# Cluster into groups
# Attrition recovery
d_mdl_attr <- d_mdl_nest %>% 
  filter(parameter %in% c("p", "r_new")) %>% 
  unnest() %>% 
  arrange(trial_id) %>% 
  mutate(
    parameter = case_match(
      parameter,
      "p" ~ "arm_risk",
      "r_new" ~ "arm_n"
    )
  ) %>% 
  rename(attrition = parameter) %>% 
  mutate(across(mean:`97.5%`, ~ round(., 2)))

# Gompertz parameters
d_mdl_gompertz <- d_mdl_nest %>% 
  filter(parameter %in% c("a", "b", "log_b", "error")) %>% 
  unnest() %>% 
  arrange(trial_id) %>% 
  mutate(
    parameter = if_else(parameter == "error", "a_error", parameter),
    parameter = factor(parameter, levels = c("log_b", "b", "a", "a_error"))
  ) %>% 
  rename(gompertz = parameter) %>% 
  arrange(trial_id, gompertz) %>% 
  mutate(across(mean:`97.5%`, ~ round(., 2)))

# Intercept
d_mdl_int <- d_mdl_nest %>% 
  filter(parameter %in% c("mu1", "mu_mu1", "tau_mu1", "z_mu1")) %>% 
  unnest() %>% 
  arrange(trial_id) %>% 
  mutate(
    parameter = factor(parameter, levels = c("mu_mu1", "tau_mu1", "mu1", "z_mu1")),
    level = case_match(
      parameter,
      "mu_mu1" ~ "global",
      "tau_mu1" ~ "global",
      .default = "local"
    ),
    level = factor(level, levels = c("global", "local"))
  ) %>% 
  arrange(level, trial_id) %>% 
  select(level, parameter, everything()) %>% 
  mutate(across(mean:`97.5%`, ~ round(., 2)))

# Beta shrinkage
d_mdl_shrink <- d_mdl_sum %>% 
  filter(parameter %in% c("lambda_beta", "tau_beta", "z_beta")) %>% 
  left_join(
    b_add_names %>% 
      distinct(tidy_name) %>% 
      filter(!is.na(tidy_name)) %>% 
      mutate(draw = seq(30))
  ) %>% 
  mutate(
    parameter = case_match(
      parameter,
      "lambda_beta" ~ "lambda",
      "tau_beta" ~ "global",
      "z_beta" ~ "latent"
    ),
    parameter = factor(parameter, levels = c("global", "latent", "lambda"))
  ) %>% 
  arrange(parameter) %>% 
  select(
    parameter, tidy_name, everything(),
    -draw
  ) %>% 
  rename(
    covariate = tidy_name
  ) %>% 
  mutate(across(mean:`97.5%`, ~ round(., 2)))

# Fit (log-likelihood, log posterior probability density)
d_mdl_loglik <- d_mdl_nest %>% 
  filter(parameter %in% c("log_lik", "lp__")) %>% 
  unnest() %>% 
  arrange(trial_id) %>% 
  mutate(
    parameter = case_match(
      parameter,
      "lp__" ~ "log-posterior density",
      "log_lik" ~ "log-likelihood"
    ),
    parameter = factor(
      parameter,
      levels = c("log-posterior density", "log-likelihood")
    )
  ) %>% 
  arrange(parameter, trial_id) %>% 
  mutate(across(mean:`97.5%`, ~ round(., 2)))

# Save
write_csv(d_mdl_attr, "outputs/manuscript_summaries/design_model/posterior_attr_recovery.csv")
write_csv(d_mdl_gompertz, "outputs/manuscript_summaries/design_model/posterior_gompertz.csv")
write_csv(d_mdl_int, "outputs/manuscript_summaries/design_model/posterior_intercept.csv")
write_csv(d_mdl_shrink, "outputs/manuscript_summaries/design_model/posterior_shrinkage.csv")
write_csv(d_mdl_loglik, "outputs/manuscript_summaries/design_model/posterior_fit_diags.csv")

# Posterior summaries: Attrition odds ratios -----------------------------------

# Create draws matrix
e_draws_mat <- as_draws_matrix(as.array(a_imp_fit))

# Extract beta columns
e_beta_cols <- grep("^beta\\[", colnames(e_draws_mat), value = TRUE)
e_beta_draws <- e_draws_mat[, e_beta_cols, drop = FALSE]

# Extract global log-rate parameter
e_mu_mu1_draws <- e_draws_mat[, "mu_mu1"]

# Median follow-up duration across trials (years)
e_dur <- median(a_imp_stan$dt, na.rm = TRUE)

# Baseline design vector (mean trial design)
e_base <- colMeans(a_imp_stan$X)
e_covars <- colnames(a_imp_stan$X)

# Identify groups of categorical variables
e_covar_groups <- list(
  e_covars[str_detect(e_covars, "^trttype")],
  e_covars[str_detect(e_covars, "^n_cat")],
  e_covars[str_detect(e_covars, "^arm_n_cat")],
  e_covars[str_detect(e_covars, "^n_sites_cat")],
  e_covars[str_detect(e_covars, "^masking")],
  e_covars[str_detect(e_covars, "^upper_bmi")],
  e_covars[str_detect(e_covars, "^lower_egfr")]
)

# Function to identify categorical contrast columns
contrast_cols <- function(var) {
  
  cols <- keep(e_covar_groups, ~ var %in% .x)
  
  if (length(cols) == 0) {
    var
  } else {
    cols[[1]]
  }
  
}

# Function to predict attrition risk
predict_risk <- function(x) {
  
  # Linear predictor
  rate_lp <- as.numeric(e_beta_draws %*% x)
  
  # Rate
  log_rate <- e_mu_mu1_draws + rate_lp
  rate <- exp(log_rate)
  
  # Shape parameter for a typical trial
  shape <- a_imp_stan$alpha_shape +
    a_imp_stan$beta_lograte * log_rate
  
  # Gompertz cumulative distribution
  risk <- pgompertz(e_dur, shape = shape, rate = rate)
  pmin(
    pmax(as.numeric(risk), .Machine$double.eps),
    1 - .Machine$double.eps
  )
}

# Estimate odds ratio for each trial characteristic
e_or <- map_dfr(seq_along(e_covars), function(i){
  
  x_new <- e_base
  x_base <- e_base
  var <- e_covars[i]
  cols <- contrast_cols(var)
  
  # Binary or categorical contrast
  x_new[cols] <- 0
  x_base[cols] <- 0
  x_new[var] <- 1
  
  risk_new <- predict_risk(x_new)
  risk_base <- predict_risk(x_base)
  odds_ratio <- exp(qlogis(risk_new) - qlogis(risk_base))
  
  tibble(
    raw_name = var,
    mean = mean(odds_ratio),
    median = median(odds_ratio),
    sd = sd(odds_ratio),
    `2.5%` = quantile(odds_ratio, .025),
    `97.5%` = quantile(odds_ratio, .975)
  )
  
})

# Add tidy covariate names and plot
e_or_tidy <- b_var_names %>% 
  distinct(raw_name, tidy_name) %>% 
  left_join(e_or)

e_or_plot <- e_or_tidy %>% 
  ggplot(aes(x = mean, xmin = `2.5%`, xmax = `97.5%`, y = fct_rev(tidy_name))) +
  geom_point() +
  geom_linerange() +
  geom_vline(xintercept = 1, colour = "red", linetype = "dashed") +
  scale_x_log10(
    limits = c(0.5, 2.5),
    breaks = c(0.5, 1, 2),
    labels = c("0.5", "1.0", "2.0"),
    minor_breaks = NULL
  ) +
  theme_classic() +
  labs(
    x = "Mean odds ratio for attrition at median follow-up (95% credible intervals)", 
    y = NULL
  )

# Save
write_csv(
  e_or_tidy %>% mutate(across(mean:`97.5%`, ~ round(., 3))), 
  "outputs/manuscript_summaries/design_model/trial_chars_odds_ratio.csv"
)

ggsave(
  "outputs/design/design_trial_chars_odds_ratio.png",
  e_or_plot,
  height = 6,
  width = 8,
  units = "in"
)

# Model validation -------------------------------------------------------------

# Import posterior draws for each fold
a_imp_stan
f_cv_imp <- readRDS("outputs/design/internal_validation/out/cv_predictions.rds")

# Function to summarise mean and credible intervals
summ_mean_q <- function(x, name, na.rm = TRUE) {
  tibble::tibble(
    "{name}"       := mean(x, na.rm = na.rm),
    "{name}_q025"  := quantile(x, 0.025, na.rm = na.rm, names = FALSE),
    "{name}_q975"  := quantile(x, 0.975, na.rm = na.rm, names = FALSE)
  )
}

# Posterior summaries
f_diag_summary <- f_diag_all %>%
  pivot_longer(mae:cor, names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  summarise(
    mean   = mean(value),
    median = median(value),
    sd     = sd(value),
    `2.5%` = quantile(value, .025),
    `97.5%` = quantile(value, .975),
    .groups = "drop"
  ) %>%
  mutate(across(mean:`97.5%`, ~round(.x, 3)))

# Calibration curve (binned)
f_calib_df <- f_cv_df %>%
  mutate(
    risk_bin = ntile(risk_pred, 10)
  ) %>%
  group_by(risk_bin) %>%
  summarise(
    pred_mean = mean(risk_pred),
    obs_mean  = mean(risk_obs),
    n = n(),
    .groups = "drop"
  )

# Calibration plot
f_cv_plot <- ggplot() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    colour = "red"
  ) +
  # raw predictions
  geom_point(
    data = f_cv_df,
    aes(x = risk_pred, y = risk_obs),
    size = 0.5,
    alpha = 0.1
  ) +
  # binned calibration curve
  geom_line(
    data = f_calib_df,
    aes(x = pred_mean, y = obs_mean),
    colour = "steelblue",
    linewidth = 0.5
  ) +
  geom_point(
    data = f_calib_df,
    aes(x = pred_mean, y = obs_mean),
    colour = "steelblue",
    size = 1
  ) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 0.9), n.breaks = 8) +
  scale_x_continuous(limits = c(0, 0.6), n.breaks = 8) +
  labs(
    x = "Predicted attrition risk",
    y = "Observed attrition risk"
  )

ggsave(
  "outputs/plots/design_valid_calibration_plot.png",
  f_cv_plot,
  width = 4,
  height = 4,
  units = "in"
)
