
# Test run
# Internal validation of design-stage model
# Uses k-fold (k = 5) cross-validation
# Script is run within MARS (HPC)

rm(list = ls())
gc()

set.seed(123)

library(tidyverse)
library(bayesplot)
library(posterior)
library(rstan)

# Rstan config
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Number of folds
K <- 5
out_dir <- "outputs/"

# Import stan dataset
full_dat <- readRDS("latest_stan_design.rds")
I_full <- full_dat$I

# Create arm-level folds
fold_id <- sample(rep(1:K, length.out = I_full))
saveRDS(fold_id, file.path(out_dir, "fold_id.rds"))

# Compile stan model
mdl <- stan_model("mdl_design.stan")

# Gompertz CDF
gompertz_prob <- function(b, a, dt) {
  ifelse(
    abs(a) < 1e-4,
    1 - exp(-b * dt),
    1 - exp(-b / a * (exp(a * dt) - 1))
  )
}

# List to store CV res
cv_results <- vector("list", K)

# K-fold loop
for (k in 1:K) {
  
  message(sprintf("Running fold %d / %d", k, K))
  
  train_idx <- which(fold_id != k)
  test_idx  <- which(fold_id == k)
  
  # Training data
  stan_dat <- full_dat
  stan_dat$I  <- length(train_idx)
  stan_dat$r  <- full_dat$r[train_idx]
  stan_dat$N  <- full_dat$N[train_idx]
  stan_dat$dt <- full_dat$dt[train_idx]
  stan_dat$lp <- full_dat$lp[train_idx]
  stan_dat$X  <- full_dat$X[train_idx, , drop = FALSE]
  
  # Fit model
  fit_k <- sampling(
    mdl,
    data = stan_dat,
    warmup = 1000,
    iter = 4000,
    chains = 4,
    seed = 123 + k,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
  
  saveRDS(fit_k, file.path(out_dir, sprintf("fit_fold_%d.rds", k)))
  
  # Extract posterior draws
  draws <- as_draws_df(fit_k)
  
  mu_mu1   <- draws$mu_mu1
  tau_mu1  <- draws$tau_mu1
  beta_mat <- as.matrix(draws %>% select(starts_with("beta[")))
  
  S <- length(mu_mu1)
  n_test <- length(test_idx)
  
  # Test-set data
  X_test  <- full_dat$X[test_idx, , drop = FALSE]
  dt_test <- full_dat$dt[test_idx]
  lp_test <- full_dat$lp[test_idx]
  N_test  <- full_dat$N[test_idx]
  r_test  <- full_dat$r[test_idx]
  
  # Posterior predictive draws (test data)
  pred_risk <- matrix(NA_real_, nrow = S, ncol = n_test)
  
  for (s in 1:S) {
    
    # Sample new random intercepts for held-out arms
    mu1_new <- rnorm(n_test, mu_mu1[s], tau_mu1[s])
    
    log_b <- mu1_new + X_test %*% beta_mat[s, ]
    b     <- exp(log_b)
    a     <- b * full_dat$beta_rate + lp_test
    
    pred_risk[s, ] <- gompertz_prob(b, a, dt_test)
  }
  
  # Summarise predictions
  pred_risk_med <- apply(pred_risk, 2, median)
  pred_count    <- N_test * pred_risk_med
  obs_risk      <- r_test / N_test
  
  cv_results[[k]] <- tibble(
    fold = k,
    arm  = test_idx,
    N    = N_test,
    r_obs = r_test,
    risk_obs  = obs_risk,
    risk_pred = pred_risk_med,
    count_pred = pred_count
  )
  
  gc()
}

# Combine cv res
cv_df <- bind_rows(cv_results)
saveRDS(cv_df, file.path(out_dir, "cv_results.rds"))

# Summarise diagnostics
cv_summary <- cv_df %>%
  summarise(
    mean_obs_risk  = mean(risk_obs),
    mean_pred_risk = mean(risk_pred),
    calib_in_large = mean(risk_pred - risk_obs),
    mae_risk       = mean(abs(risk_pred - risk_obs))
  )

write_csv(cv_summary, file.path(out_dir, "cv_summary.csv"))

# Calibration plot
cv_plot <- ggplot(cv_df, aes(x = risk_obs, y = risk_pred)) +
  geom_abline(linetype = 2) +
  geom_point(alpha = 0.6) +
  labs(
    x = "Observed attrition risk",
    y = "Predicted attrition risk"
  ) +
  theme_classic()

ggsave(
  file.path(out_dir, "cv_calibration_plot.png"),
  cv_plot,
  width = 6,
  height = 6,
  units = "in"
)
