
# Interal validation
# Held-out cross-validation to evaluate predictive performance for trials
# with >1 arm
# HOCV model fit with trials with 1 arm
# Posterior draws used to predict attrition probability in held-out trials
# Evaluate predictive performance using Brier score
# Arm-level calibration plot and Brier score for full model after HOCV

source("scripts/00_config.R")
library(tidyverse)
library(bayesplot)
library(tidybayes)
library(posterior)
library(rstan)
library(loo)

# Import fitted model
a_imp_df <- readRDS("processed_data/latest_dataset.rds")
a_imp_mdl <- readRDS("mars/design_model/outputs/design_fit.rds")
a_imp_draws <- readRDS("mars/design_model/outputs/design_draws.rds")

# Get draws df and save for future use (imported above)
# b_draws_df <- as_draws_df(a_imp_mdl)
# saveRDS(b_draws_df, "mars/design_model/outputs/design_draws.rds")
# b_draws_df_hocv <- as_draws_df(a_imp_mdl_hocv)
# saveRDS(b_draws_df_hocv, "mars/design_hocv/outputs/design_hocv_draws.rds")

# LOO for HOCV model - Fails, most are k > 0.7
# Issues caused by hierarchical structure, highly influential observations
# will cause it to fail if removed, dependence on >1 hierarchical prior
# a_loo_loglik <- extract_log_lik(a_imp_mdl_hocv, "log_lik", merge_chains = FALSE)
# a_loo_eff <- relative_eff(exp(a_loo_loglik), cores = 2)
# a_loo <- loo(a_loo_loglik, r_eff = a_loo_eff, cores = 2, method = "sis")
# print(a_loo)
# plot(a_loo)

# Hold-out cross-validation ----------------------------------------------------

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 1. Get held-out data and posterior draw parameters
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Held-out covariate matrix, arm_n, duration and attrition counts
b_ho_mdl_mat <- a_imp_df_hocv$X
b_ho_arm_n <- a_imp_df_hocv$N
b_ho_dur <- a_imp_df_hocv$dt
b_ho_attr <- a_imp_df_hocv$r

# Create matrix to store predicted probabilities
b_ho_sample <- 12000
b_ho_n_arms <- a_imp_df_hocv$I
b_ho_pred_attr <- matrix(NA, nrow = b_ho_n_arms, ncol = b_ho_sample)

# Extract relevant parameters from posterior draws
b_param_draws <- a_imp_draws_hocv %>% 
  pivot_longer(everything()) %>% 
  filter(grepl("mu|tau|beta", name))

# Tidy parameters
b_param_tidy <- b_param_draws %>% 
  separate(name, into = c("param", "draw"), sep = "\\[") %>% 
  mutate(row = gsub("\\]", "", draw))

# 3,996,000 = 366 * 12000 (arm-level parameters)
# 228,000 = 19 * 12000 (coefficients and latent variables)
# 12000 (global means and sds for mu1 and beta)
b_param_tidy %>% count(param)

# Split by parameter type
b_param_split <- split(b_param_tidy, b_param_tidy$param)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 2. Compute mu1
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Sample new z_mu1 values for held-out arms
set.seed(123)
z_mu1_new <- matrix(rnorm(12000*33, 0, 1), nrow = 12000, ncol = 33)

# Compute mu1 for held-out arms
b_ho_mu1 <- b_param_split$mu_mu1$value + b_param_split$tau_mu1$value * z_mu1_new
dim(b_ho_mu1) # 12000x33


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 3. Compute betas
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Wide format for beta parameter draws
b_ho_beta_params <- b_param_split$z_beta %>% 
  mutate(row = as.integer(row)) %>% 
  select(-draw) %>% 
  pivot_wider(names_from = param, values_from = row)

b_ho_beta_lambda_params <- b_param_split$lambda_beta %>% 
  mutate(row = as.integer(row)) %>% 
  pivot_wider(names_from = param, values_from = row) %>%
  select(-draw)

b_ho_beta_tau <- b_param_split$tau_beta$value

# Construct matrices for beta parameter draws
z_beta <- matrix(b_ho_beta_params$value, nrow = 12000, ncol = 19)
beta_lambda <- matrix(b_ho_beta_lambda_params$value, nrow = 12000, ncol = 19)
b_ho_beta <- z_beta * beta_lambda * b_ho_beta_tau
dim(b_ho_beta) # 12000x19


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 4. Compute rate and shape
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Compute linear predictor
b_ho_rate_lp <- b_ho_mdl_mat %*% t(b_ho_beta)
b_ho_rate_lp <- t(b_ho_rate_lp)

# Compute rate
b_ho_rate_mat <- exp(b_ho_mu1 + b_ho_rate_lp)
dim(b_ho_rate_mat) # 12000x33

# Compute shape using rate plus product of covariates * betas
b_ho_shape_beta <- a_imp_stan_hocv$B
b_ho_shape_mat <- b_ho_rate_mat + matrix(
  b_ho_mdl_mat %*% b_ho_shape_beta, 
  nrow = 12000, ncol = 33, byrow = TRUE
)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 4. Compute predicted attrition using Gompertz CDF
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Matrix for duration
b_ho_dur_mat <- matrix(rep(b_ho_dur, each = b_ho_sample), nrow = b_ho_sample)

# Set conditions for using exponential
eps <- 1e-4
is_exp <- abs(b_ho_shape_mat) < eps

# Exponential cases
b_ho_pred_attr[is_exp] <- 1 - exp(-b_ho_rate_mat[is_exp] * b_ho_dur_mat[is_exp])

# Gompertz cases
is_gomp <- which(!is_exp)
b_ho_pred_attr[is_gomp] <- 1 - exp(
  -b_ho_rate_mat[is_gomp] / b_ho_shape_mat[is_gomp] *
    (exp(b_ho_shape_mat[is_gomp] * b_ho_dur_mat[is_gomp]) - 1)
)

dim(b_ho_pred_attr)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 5. Get attrition predictions and tidy summary
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Get tidy dataset of held-out attrition predictions with groups
b_ho_tidy <- as_tibble(t(b_ho_pred_attr)) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(V1:V33)
  
# Get names
b_ho_tidy_names <- a_imp_ho_ids %>% 
  mutate(name = paste0("V", seq(1:33))) %>% 
  group_by(trial_id, class, trttype) %>% 
  mutate(row = row_number()) %>% 
  ungroup()

# Add names
b_ho_names <- b_ho_tidy_names %>% left_join(b_ho_tidy)

# Get summary
b_ho_sum <- b_ho_names %>% 
  group_by(trial_id, trttype, class, row) %>% 
  reframe(
    mean = mean(value),
    q2.5 = quantile(value, 0.025), 
    q97.5 = quantile(value, 0.975)
  )

# Add observed attrition
b_ho_obs <- b_ho_sum %>% 
  left_join(
    a_imp_df %>% 
      select(trial_id, class, trttype, arm_n, attr, attr_prop) %>% 
      group_by(trial_id, class, trttype) %>% 
      mutate(row = row_number()) %>% 
      ungroup()
  )

# Get summary values as matrix
b_ho_sum <- apply(
  b_ho_pred_attr, 2,
  \(x) c(mean = mean(x), q2.5 = quantile(x, 0.025), q97.5 = quantile(x, 0.975))
)

# Get tidy summary
b_ho_sum_tidy <- as_tibble(t(b_ho_sum))

# Brier score ------------------------------------------------------------------

# Get predicted attrition
b_probs <- a_imp_draws %>% 
  pivot_longer(everything()) %>% 
  filter(grepl("^p", name)) %>% 
  separate(name, into = c("var", "row"), sep = "\\[") %>% 
  mutate(
    row = gsub("\\]", "", row),
    row = as.integer(row),
    iter = rep(1:(n() / 366), each = 366)
  )

b_id <- a_imp_df %>% 
  select(trial_id, trttype, arm_n, attr, attr_prop) %>% 
  mutate(row = as.integer(row_number()))

# Join observed variables
b_join <- b_id %>%
  left_join(b_probs) %>% 
  rename(
    obs_attr = attr_prop,
    pred_attr = value
  ) %>% 
  select(-var) %>% 
  mutate(
    trttype = case_when(
      trttype == 0L ~ "Plc",
      TRUE ~ "AC"
    ),
    trttype = factor(trttype, levels = c("Plc", "AC"))
  )

# Compute Brier score
b_brier <- b_join %>% 
  group_by(trttype, iter) %>% 
  reframe(brier = mean((obs_attr - pred_attr)^2))

# Summarise Brier scores across posterior draws
b_brier_sum <- b_brier %>% 
  group_by(trttype) %>% 
  reframe(
    mean_brier = mean(brier),
    sd_brier = sd(brier),
    med_brier = median(brier),
    min_brier = min(brier),
    max_brier = max(brier)
  )

write_csv(b_brier_sum, "mars/design_model/outputs/brier_sum.csv")

# Plot Brief scores
b_brier_dist <- b_brier %>%
  mutate(trttype = as.factor(trttype)) %>% 
  ggplot(aes(x = brier, group = trttype)) +
  geom_density(aes(fill = trttype), alpha = 0.3, colour = "black") +
  scale_x_continuous(n.breaks = 8) +
  theme_classic() +
  labs(x = "Brier score", y = "Density", fill = "Treatment type")

b_brier_dist

ggsave(
  "mars/design_model/outputs/brier_dist.png",
  b_brier_dist,
  width = 8,
  height = 8,
  units = "in"
)

# Calibration plot -------------------------------------------------------------

# Summarise predicted attrition
# Create bins of predicted attrition
c_calib_df <- b_join %>% 
  group_by(row, trttype) %>% 
  reframe(
    obs_attr = first(obs_attr),
    pred_attr = mean(pred_attr)
  ) %>% 
  mutate(
    bin = ntile(pred_attr, 10)
  ) %>% 
  group_by(bin, trttype) %>% 
  reframe(
    mean_pred = mean(pred_attr),
    mean_obs = mean(obs_attr),
    n = n()
  )

# Plot
c_calib_plot <- c_calib_df %>% 
  ggplot(aes(x = mean_pred, y = mean_obs, group = trttype)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  geom_point() +
  geom_line() +
  facet_wrap(~ trttype, scales = "free") +
  theme_classic() +
  labs(x = "Mean predicted attrition", y = "Observed attrition")

c_calib_plot

ggsave(
  "mars/design_model/outputs/calib_plot.png",
  c_calib_plot,
  width = 8,
  height = 8,
  units = "in"
)