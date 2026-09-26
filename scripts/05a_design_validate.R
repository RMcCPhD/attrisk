
# Internally validate design-stage model
# K-fold cross-validation performed
# Script is run on a high performance computing unit
# Can run locally but takes an extremely long time due to computing demands

rm(list = ls())
gc()
set.seed(123)

# Packages
library(tidyverse)
library(rstan)
library(bayesplot)
library(posterior)

# Config
options(scipen = 999)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Define folds -----------------------------------------------------------------

# Set number of folds and save directory
K_folds <- 10
out_dir <- "kfold_cv"

# Import Stan data
imp_dat <- readRDS("stan_data.rds")

# Create folds
# Folds are trial-level since arms are collapsed where >1 non-experimental
# Stratified by observed risk to balance low and high-risk trials, i.e.
# minimise chance of all extremes getting bundled together
fold_tbl <- tibble(
  row = seq_len(imp_dat$I),
  r = imp_dat$r,
  N = imp_dat$N,
  risk_obs = r / N,
  fold = sample(rep(seq_len(K_folds), length.out = imp_dat$I))
)

folds <- fold_tbl$fold

# Save fold IDs
saveRDS(folds, file.path(out_dir, "fold_id.rds"))
write_csv(fold_tbl, file.path(out_dir, "fold_assignment.csv"))

# Summarise trials within folds
folds_sum <- fold_tbl %>% 
  group_by(fold) %>% 
  reframe(
    n = n(),
    mean_obs_risk = mean(risk_obs),
    median_obs_risk = median(risk_obs),
    min_obs_risk = min(risk_obs),
    max_obs_risk = max(risk_obs),
    total_r = sum(r),
    total_N = sum(N),
    aggregate_obs_risk = total_r / total_N
  )

# Save summary
write_csv(folds_sum, file.path(out_dir, "fold_summary.csv"))

# Compile model
compile_mdl <- stan_model("design_model.stan")

# Helper functions -------------------------------------------------------------

# Gompertz CDF
# Numerical stability baked in, reduces to exponential if shape is near-zero
gompertz_prob <- function(b, a, dt) {
  ifelse(
    abs(a) < 1e-4,
    1 - exp(-b * dt),
    1 - exp(-b / a * (exp(a * dt) - 1))
  )
}

# Stable log-mean-exp for held-out log predictive density
log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

# Subset Stan data by row-indexed components only
# Avoids having to manually update every element and protects global parameters
subset_stan_rows <- function(dat, idx) {
  
  dat_out <- dat
  
  old_I <- dat$I
  new_I <- length(idx)
  
  dat_out <- purrr::map(dat_out, function(x) {
    
    if (is.matrix(x) && nrow(x) == old_I) {
      x[idx, , drop = FALSE]
      
    } else if (is.atomic(x) && length(x) == old_I) {
      x[idx]
      
    } else {
      x
    }
  })
  
  dat_out$I <- new_I
  
  # Ensure count variables remain integer arrays for Stan
  dat_out$r <- as.integer(dat_out$r)
  dat_out$N <- as.integer(dat_out$N)
  
  dat_out
}

# Extract simple diagnostics for each fold after model fit
get_hmc_diagnostics <- function(fit, fold, train_idx, test_idx) {
  
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  
  n_divergent <- sum(
    purrr::map_int(
      sampler_params,
      ~ sum(.x[, "divergent__"])
    )
  )
  
  n_max_treedepth <- sum(
    purrr::map_int(
      sampler_params,
      ~ sum(.x[, "treedepth__"] >= 15)
    )
  )
  
  fit_sum <- summary(fit)$summary
  
  tibble(
    fold = fold,
    n_train = length(train_idx),
    n_test = length(test_idx),
    n_divergent = n_divergent,
    n_max_treedepth = n_max_treedepth,
    max_rhat = max(fit_sum[, "Rhat"], na.rm = TRUE),
    min_n_eff = min(fit_sum[, "n_eff"], na.rm = TRUE)
  )
  
}

# Cross-validation -------------------------------------------------------------

cv_diag <- vector("list", K_folds)

for (k in seq_len(K_folds)) {
  
  message("------------------------------------------------------------")
  message(sprintf("Running fold %d / %d", k, K_folds))
  message("------------------------------------------------------------")
  
  train_idx <- which(folds != k)
  test_idx  <- which(folds == k)
  
  # Save indices so local code can reconstruct held-out sets exactly
  saveRDS(
    list(
      fold = k,
      train_idx = train_idx,
      test_idx = test_idx
    ),
    file.path(out_dir, sprintf("fold_%d_indices.rds", k))
  )
  
  # Prepare training Stan data
  stan_dat_k <- subset_stan_rows(imp_dat, train_idx)
  
  # Save training data for transparency/reproducibility
  saveRDS(
    stan_dat_k,
    file.path(out_dir, sprintf("stan_data_fold_%d.rds", k))
  )
  
  # Fit model
  fit_k <- sampling(
    compile_mdl,
    data = stan_dat_k,
    warmup = 1000,
    iter = 4000,
    chains = 4,
    seed = 123 + k,
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 15
    )
  )
  
  # Save fitted model for local machine
  saveRDS(
    fit_k,
    file.path(out_dir, sprintf("fit_fold_%d.rds", k))
  )
  
  # Save fold diagnostics
  cv_diag[[k]] <- get_hmc_diagnostics(
    fit = fit_k,
    fold = k,
    train_idx = train_idx,
    test_idx = test_idx
  )
  
  write_csv(
    cv_diag[[k]],
    file.path(out_dir, sprintf("diagnostics_fold_%d.csv", k))
  )
  
  # Save cumulative diagnostics after each fold, so progress is preserved
  write_csv(
    bind_rows(cv_diag),
    file.path(out_dir, "cv_diagnostics_running.csv")
  )
  
  # Clear large objects before next fold
  rm(fit_k, stan_dat_k)
  gc()
}

# Save final diagnostics --------------------------------------------------------

# Combine fold diagnostics and save
cv_diag_df <- bind_rows(cv_diag)

write_csv(
  cv_diag_df,
  file.path(out_dir, "cv_diagnostics.csv")
)

saveRDS(
  cv_diag_df,
  file.path(out_dir, "cv_diagnostics.rds")
)