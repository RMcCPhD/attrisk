
# Generate outputs from internal cross-validation
# Summaries and plots

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

set.seed(123)

# Define file directory for folds
cv_dir <- "outputs/design/internal_validation/models/"

# Import Stan data
imp_stan <- readRDS("processed_data/stan_data.rds")

# Import fold metadata
imp_fold_ids <- readRDS("outputs/design/internal_validation/meta/fold_id.rds")

# Infer K from the saved allocation so all fitted folds are processed
fold_values <- sort(unique(imp_fold_ids))

if (length(imp_fold_ids) != imp_stan$I) {
  stop("Fold IDs do not match the number of trials in the Stan data.")
}

if (anyNA(fold_values) ||
    !setequal(fold_values, seq_len(length(fold_values)))) {
  stop("Fold IDs must be contiguous integers starting at 1.")
}

K_folds <- length(fold_values)

# Import cv diagnostics
imp_cv_diag <- read_csv("outputs/design/internal_validation/sum/cv_diagnostics.csv")

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

# Helper functions to summarise fold draws
q025 <- function(x) {
  quantile(x, 0.025, na.rm = TRUE, names = FALSE)
}

q975 <- function(x) {
  quantile(x, 0.975, na.rm = TRUE, names = FALSE)
}

col_q025 <- function(x) {
  apply(x, 2, q025)
}

col_q975 <- function(x) {
  apply(x, 2, q975)
}

safe_cor <- function(x, y) {
  ok <- complete.cases(x, y)
  
  if (sum(ok) < 2 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(NA_real_)
  }
  
  cor(x[ok], y[ok])
}

make_draw_df <- function(
    prediction,
    fold,
    test_idx,
    N_test,
    r_test,
    risk_mat,
    count_exp_mat,
    count_new_mat,
    log_lik_mat
) {
  
  S <- nrow(risk_mat)
  n_test <- ncol(risk_mat)
  
  tibble(
    prediction = prediction,
    fold = fold,
    draw = rep(seq_len(S), times = n_test),
    row = rep(test_idx, each = S),
    N = rep(N_test, each = S),
    r_obs = rep(r_test, each = S),
    risk_obs = rep(r_test / N_test, each = S),
    risk_pred = as.vector(risk_mat),
    count_pred = as.vector(count_exp_mat),
    count_new = as.vector(count_new_mat),
    log_lik = as.vector(log_lik_mat)
  )
}

# Extract fold predictions
get_fold_pred <- function(k, imp_dat, folds, cv_dir, max_draws = Inf) {
  
  # Print message indicating progress
  message(sprintf("Processing fold %d", k))
  
  # Identify held-out trials
  test_idx <- which(folds == k)

  if (length(test_idx) == 0) {
    stop(sprintf("Fold %d does not contain any held-out trials.", k))
  }
  
  # Load fold model
  fit_path <- file.path(cv_dir, sprintf("fit_fold_%d.rds", k))

  if (!file.exists(fit_path)) {
    stop("Could not find fitted fold model: ", fit_path)
  }

  fit_k <- readRDS(fit_path)
  
  # Extract prediction parameters
  ext <- extract(fit_k, pars = c("mu_mu1", "tau_mu1", "beta"))
  S_all <- length(ext$mu_mu1)
  
  # Use all posterior draws unless max_draws is specified
  if (is.finite(max_draws) && max_draws < S_all) {
    draw_idx <- sample(seq_len(S_all), size = max_draws)
  } else {
    draw_idx <- seq_len(S_all)
  }
  
  # Define parameters
  mu_mu1 <- ext$mu_mu1[draw_idx]
  tau_mu1 <- ext$tau_mu1[draw_idx]
  beta_arr <- ext$beta[draw_idx, , drop = FALSE]
  
  S <- length(draw_idx)
  n_test <- length(test_idx)
  
  # Held-out test data
  X_test <- imp_dat$X[test_idx, , drop = FALSE]
  dt_test <- imp_dat$dt[test_idx]
  N_test <- imp_dat$N[test_idx]
  r_test <- imp_dat$r[test_idx]
  
  # Shape exploratory analysis parameters
  alpha_shape <- imp_dat$alpha_shape
  beta_shape_lograte <- imp_dat$beta_lograte
  resid_sd <- imp_dat$resid_sd
  
  # Store marginal predictions (a new trial drawn from same trial population)
  pred_risk_marginal <- matrix(NA_real_, nrow = S, ncol = n_test)
  pred_count_exp_marginal <- matrix(NA_real_, nrow = S, ncol = n_test)
  pred_count_new_marginal <- matrix(NA_real_, nrow = S, ncol = n_test)
  log_lik_marginal <- matrix(NA_real_, nrow = S, ncol = n_test)
  
  # Store typical predictions (prediction deviations set to central values)
  pred_risk_typical <- matrix(NA_real_, nrow = S, ncol = n_test)
  pred_count_exp_typical <- matrix(NA_real_, nrow = S, ncol = n_test)
  pred_count_new_typical <- matrix(NA_real_, nrow = S, ncol = n_test)
  log_lik_typical <- matrix(NA_real_, nrow = S, ncol = n_test)
  
  # Posterior marginal prediction loop
  for (s in seq_len(S)) {
    
    # Marginal predictions
    mu1_new <- rnorm(
      n = n_test,
      mean = mu_mu1[s],
      sd = tau_mu1[s]
    )
    
    log_b <- as.vector(mu1_new + X_test %*% beta_arr[s, ])
    b <- exp(log_b)
    
    z_shape_new <- rnorm(n_test, mean = 0, sd = 1)
    a <- alpha_shape + beta_shape_lograte * log_b + resid_sd * z_shape_new
    
    p <- gompertz_prob(b = b, a = a, dt = dt_test)
    
    pred_risk_marginal[s, ] <- p
    pred_count_exp_marginal[s, ] <- N_test * p
    pred_count_new_marginal[s, ] <- rbinom(n_test, size = N_test, prob = p)
    log_lik_marginal[s, ] <- dbinom(r_test, size = N_test, prob = p, log = TRUE)
    
    # Typical predictions
    mu1_typical <- rep(mu_mu1[s], n_test)
    
    log_b_typical <- as.vector(mu1_typical + X_test %*% beta_arr[s, ])
    b_typical <- exp(log_b_typical)
    
    a_typical <- alpha_shape + beta_shape_lograte * log_b_typical
    
    p_typical <- gompertz_prob(
      b = b_typical,
      a = a_typical,
      dt = dt_test
    )
    
    pred_risk_typical[s, ] <- p_typical
    pred_count_exp_typical[s, ] <- N_test * p_typical
    pred_count_new_typical[s, ] <- rbinom(n_test, size = N_test, prob = p_typical)
    log_lik_typical[s, ] <- dbinom(
      r_test,
      size = N_test,
      prob = p_typical,
      log = TRUE
    )
    
  }
  
  # Convert posterior count draws to risk
  pred_risk_new_marginal <- sweep(pred_count_new_marginal, 2, N_test, "/")
  pred_risk_new_typical <- sweep(pred_count_new_typical, 2, N_test, "/")
  
  # Row-level held-out prediction summaries
  pred_sum <- tibble(
    fold = k,
    row = test_idx,
    N = N_test,
    r_obs = r_test,
    risk_obs = r_test / N_test,
    dt = dt_test,
    
    # Marginal risk: latent attrition probability p
    risk_pred_mean = colMeans(pred_risk_marginal),
    risk_pred_median = apply(pred_risk_marginal, 2, median),
    risk_pred_q025 = col_q025(pred_risk_marginal),
    risk_pred_q975 = col_q975(pred_risk_marginal),
    
    # Marginal observed-risk posterior predictive interval r_new / N
    risk_pp_q025 = col_q025(pred_risk_new_marginal),
    risk_pp_q975 = col_q975(pred_risk_new_marginal),
    
    # Marginal expected counts
    count_pred_mean = colMeans(pred_count_exp_marginal),
    count_pred_median = apply(pred_count_exp_marginal, 2, median),
    
    # Marginal posterior predictive count intervals
    count_pp_q025 = col_q025(pred_count_new_marginal),
    count_pp_q975 = col_q975(pred_count_new_marginal),
    
    # Marginal held-out log predictive density
    elpd_i = apply(log_lik_marginal, 2, log_mean_exp),
    
    # Typical risk: latent attrition probability p
    risk_typical_mean = colMeans(pred_risk_typical),
    risk_typical_median = apply(pred_risk_typical, 2, median),
    risk_typical_q025 = col_q025(pred_risk_typical),
    risk_typical_q975 = col_q975(pred_risk_typical),
    
    # Typical observed-risk posterior predictive interval r_new / N
    risk_typical_pp_q025 = col_q025(pred_risk_new_typical),
    risk_typical_pp_q975 = col_q975(pred_risk_new_typical),
    
    # Typical expected counts
    count_typical_mean = colMeans(pred_count_exp_typical),
    count_typical_median = apply(pred_count_exp_typical, 2, median),
    
    # Typical posterior predictive count intervals
    count_typical_pp_q025 = col_q025(pred_count_new_typical),
    count_typical_pp_q975 = col_q975(pred_count_new_typical),
    
    # Typical held-out log predictive density
    elpd_i_typical = apply(log_lik_typical, 2, log_mean_exp)
  )
  
  # Draw-level predictions for credible intervals around summary metrics
  pred_draws <- bind_rows(
    make_draw_df(
      prediction = "marginal",
      fold = k,
      test_idx = test_idx,
      N_test = N_test,
      r_test = r_test,
      risk_mat = pred_risk_marginal,
      count_exp_mat = pred_count_exp_marginal,
      count_new_mat = pred_count_new_marginal,
      log_lik_mat = log_lik_marginal
    ),
    make_draw_df(
      prediction = "typical",
      fold = k,
      test_idx = test_idx,
      N_test = N_test,
      r_test = r_test,
      risk_mat = pred_risk_typical,
      count_exp_mat = pred_count_exp_typical,
      count_new_mat = pred_count_new_typical,
      log_lik_mat = log_lik_typical
    )
  )
  
  list(
    pred = pred_sum,
    draws = pred_draws
  )
  
}

# Generate cross-validation predictions ----------------------------------------

cv_out <- map(
  seq_len(K_folds),
  ~ get_fold_pred(
    k = .x,
    imp_dat = imp_stan,
    folds = imp_fold_ids,
    cv_dir = cv_dir,
    max_draws = Inf
  )
)

cv_df <- map_dfr(cv_out, "pred")
cv_draws <- map_dfr(cv_out, "draws")

# Save
saveRDS(cv_df, "outputs/design/internal_validation/out/cv_predictions.rds")
saveRDS(cv_draws, "outputs/design/internal_validation/out/cv_prediction_draws.rds")

# Performance summaries --------------------------------------------------------

# Functions to summarise performance for marginal and typical predictions
# Marginal: New trial
# Typical: Predictions without new-trial random-effect uncertainty

posterior_metric_summary <- function(x) {
  tibble(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    q975 = quantile(x, 0.975, na.rm = TRUE, names = FALSE)
  )
}

rename_summary_cols <- function(dat, prediction = c("marginal", "typical")) {
  
  prediction <- match.arg(prediction)
  
  dat %>%
    rename_with(
      .fn = ~ paste0(prediction, "_", .x),
      .cols = c(mean, median, sd, q025, q975)
    )
}

make_point_only_rows <- function(dat, prediction = c("marginal", "typical")) {
  
  prediction <- match.arg(prediction)
  
  vars <- switch(
    prediction,
    marginal = list(
      risk_pp_q025 = "risk_pp_q025",
      risk_pp_q975 = "risk_pp_q975",
      count_pp_q025 = "count_pp_q025",
      count_pp_q975 = "count_pp_q975",
      elpd = "elpd_i"
    ),
    typical = list(
      risk_pp_q025 = "risk_typical_pp_q025",
      risk_pp_q975 = "risk_typical_pp_q975",
      count_pp_q025 = "count_typical_pp_q025",
      count_pp_q975 = "count_typical_pp_q975",
      elpd = "elpd_i_typical"
    )
  )
  
  dat %>%
    transmute(
      N,
      r_obs,
      risk_obs,
      risk_pp_q025 = .data[[vars$risk_pp_q025]],
      risk_pp_q975 = .data[[vars$risk_pp_q975]],
      count_pp_q025 = .data[[vars$count_pp_q025]],
      count_pp_q975 = .data[[vars$count_pp_q975]],
      elpd_i = .data[[vars$elpd]]
    ) %>%
    mutate(
      count_covered_95 = r_obs >= count_pp_q025 & r_obs <= count_pp_q975,
      risk_pp_covered_95 = risk_obs >= risk_pp_q025 & risk_obs <= risk_pp_q975
    ) %>%
    reframe(
      total_elpd = sum(elpd_i),
      coverage_95_risk_pp = mean(risk_pp_covered_95),
      coverage_95_count = mean(count_covered_95)
    ) %>%
    pivot_longer(everything(), names_to = "sum", values_to = "mean") %>%
    mutate(
      median = NA_real_,
      sd = NA_real_,
      q025 = NA_real_,
      q975 = NA_real_
    ) %>%
    rename_summary_cols(prediction = prediction)
}

make_draw_summary <- function(draws, prediction = c("marginal", "typical")) {
  
  prediction <- match.arg(prediction)
  
  draws %>%
    filter(.data$prediction == .env$prediction) %>%
    mutate(
      # Count-scale expected prediction error
      # Use E(r | training data) rather than a random posterior predictive
      # replicate, which adds binomial sampling error to the comparison
      count_error = count_pred - r_obs,
      abs_count_error = abs(count_error),
      sq_count_error = count_error^2,
      rel_count_error = count_error / pmax(r_obs, 1),
      abs_rel_count_error = abs(rel_count_error),
      
      # Risk-scale expected prediction error.
      # These metrics compare latent predicted attrition probabilities against
      # observed attrition risks and are retained for the existing CV summaries.
      risk_error = risk_pred - risk_obs,
      abs_risk_error = abs(risk_error),
      sq_risk_error = risk_error^2,
      rel_risk_error = risk_error / pmax(risk_obs, 1e-6),
      abs_rel_risk_error = abs(rel_risk_error)
    ) %>%
    group_by(draw) %>%
    reframe(
      n = n(),
      
      mean_obs_risk = mean(risk_obs),
      mean_pred_risk = mean(risk_pred),
      mean_calibration = mean(risk_error),
      
      aggregate_obs_risk = sum(r_obs) / sum(N),
      aggregate_pred_risk = sum(count_pred) / sum(N),
      aggregate_calibration = sum(count_pred) / sum(N) - sum(r_obs) / sum(N),
      
      # Risk-scale predictive accuracy and association
      mae_risk = mean(abs_risk_error),
      rmse_risk = sqrt(mean(sq_risk_error)),
      mean_rel_error_risk = mean(rel_risk_error),
      mean_abs_rel_error_risk = mean(abs_rel_risk_error),
      median_abs_rel_error_risk = median(abs_rel_risk_error),
      cor_risk = safe_cor(risk_obs, risk_pred),
      
      # Count-scale posterior predictive recovery
      mpe_count = mean(count_error),
      mae_count = mean(abs_count_error),
      rmse_count = sqrt(mean(sq_count_error)),
      mean_rel_error_count = mean(rel_count_error),
      mean_abs_rel_error_count = mean(abs_rel_count_error),
      median_abs_rel_error_count = median(abs_rel_count_error),
      cor_count = safe_cor(r_obs, count_pred)
    ) %>%
    select(-draw) %>%
    pivot_longer(everything(), names_to = "sum", values_to = "value") %>%
    group_by(sum) %>%
    reframe(posterior_metric_summary(value)) %>%
    rename_summary_cols(prediction = prediction)
}

cv_sum_marginal <- bind_rows(
  make_draw_summary(cv_draws, prediction = "marginal"),
  make_point_only_rows(cv_df, prediction = "marginal")
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cv_sum_typical <- bind_rows(
  make_draw_summary(cv_draws, prediction = "typical"),
  make_point_only_rows(cv_df, prediction = "typical")
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cv_sum <- cv_sum_marginal %>%
  left_join(cv_sum_typical, by = "sum") %>%
  select(
    sum,
    marginal_mean,
    marginal_median,
    marginal_sd,
    marginal_q025,
    marginal_q975,
    typical_mean,
    typical_median,
    typical_sd,
    typical_q025,
    typical_q975
  )

write_csv(
  cv_sum,
  "outputs/design/internal_validation/out/cv_sum.csv"
)

# Fold-specific summaries ------------------------------------------------------

make_fold_draw_summary <- function(draws, prediction = c("marginal", "typical")) {
  
  prediction <- match.arg(prediction)
  
  draws %>%
    filter(.data$prediction == .env$prediction) %>%
    mutate(
      risk_error = risk_pred - risk_obs,
      abs_risk_error = abs(risk_error),
      sq_risk_error = risk_error^2,
      
      # Count-scale expected prediction error
      count_error = count_pred - r_obs,
      abs_count_error = abs(count_error),
      sq_count_error = count_error^2
    ) %>%
    group_by(fold, draw) %>%
    reframe(
      n = n(),
      mean_obs_risk = mean(risk_obs),
      mean_pred_risk = mean(risk_pred),
      aggregate_obs_risk = sum(r_obs) / sum(N),
      aggregate_pred_risk = sum(count_pred) / sum(N),
      aggregate_calibration = sum(count_pred) / sum(N) - sum(r_obs) / sum(N),
      mae_risk = mean(abs_risk_error),
      rmse_risk = sqrt(mean(sq_risk_error)),
      mpe_count = mean(count_error),
      mae_count = mean(abs_count_error),
      rmse_count = sqrt(mean(sq_count_error)),
      cor_count = safe_cor(r_obs, count_pred)
    ) %>%
    select(-draw) %>%
    pivot_longer(
      cols = -fold,
      names_to = "sum",
      values_to = "value"
    ) %>%
    group_by(fold, sum) %>%
    reframe(posterior_metric_summary(value)) %>%
    rename_summary_cols(prediction = prediction)
}

make_fold_point_only_rows <- function(dat, prediction = c("marginal", "typical")) {
  
  prediction <- match.arg(prediction)
  
  vars <- switch(
    prediction,
    marginal = list(elpd = "elpd_i"),
    typical = list(elpd = "elpd_i_typical")
  )
  
  dat %>%
    group_by(fold) %>%
    reframe(total_elpd = sum(.data[[vars$elpd]])) %>%
    pivot_longer(
      cols = -fold,
      names_to = "sum",
      values_to = "mean"
    ) %>%
    mutate(
      median = NA_real_,
      sd = NA_real_,
      q025 = NA_real_,
      q975 = NA_real_
    ) %>%
    rename_summary_cols(prediction = prediction)
}

cv_sum_fold_marginal <- bind_rows(
  make_fold_draw_summary(cv_draws, prediction = "marginal"),
  make_fold_point_only_rows(cv_df, prediction = "marginal")
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cv_sum_fold_typical <- bind_rows(
  make_fold_draw_summary(cv_draws, prediction = "typical"),
  make_fold_point_only_rows(cv_df, prediction = "typical")
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cv_sum_fold <- cv_sum_fold_marginal %>%
  left_join(cv_sum_fold_typical, by = c("fold", "sum")) %>%
  select(
    fold,
    sum,
    marginal_mean,
    marginal_median,
    marginal_sd,
    marginal_q025,
    marginal_q975,
    typical_mean,
    typical_median,
    typical_sd,
    typical_q025,
    typical_q975
  )

write_csv(
  cv_sum_fold,
  "outputs/design/internal_validation/out/cv_sum_folds.csv"
)

# Calibration by predicted risk ------------------------------------------------

# Bin observed risk in AgD
cv_calib_bins <- cv_df %>%
  mutate(pred_bin = ntile(risk_pred_mean, 10)) %>%
  group_by(pred_bin) %>%
  reframe(
    n_trials = n(),
    mean_pred_risk = mean(risk_pred_mean),
    aggregate_pred_risk = sum(count_pred_mean) / sum(N),
    mean_obs_risk = mean(risk_obs),
    aggregate_obs_risk = sum(r_obs) / sum(N),
    total_r = sum(r_obs),
    total_N = sum(N)
  ) %>% 
  mutate(across(everything(), ~ round(., 3)))

write_csv(
  cv_calib_bins,
  "outputs/design/internal_validation/out/cv_calibration_binned_risk.csv"
)

# Calibration plot
cv_calib_plot <- cv_calib_bins %>%
  ggplot(aes(x = aggregate_pred_risk, y = aggregate_obs_risk)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 2) +
  geom_line() +
  theme_classic() +
  labs(
    x = "Predicted attrition risk",
    y = "Observed attrition risk"
  )

ggsave(
  "outputs/design/internal_validation/out/cv_calibration_plot.png",
  cv_calib_plot,
  width = 6,
  height = 5,
  units = "in"
)

# Observed versus predicted risk -----------------------------------------------

cv_risk_plot <- cv_df %>%
  ggplot(aes(x = risk_obs, y = risk_pred_mean)) +
  geom_abline(linetype = "dashed") +
  geom_errorbar(
    aes(ymin = risk_pp_q025, ymax = risk_pp_q975),
    alpha = 0.25,
    width = 0
  ) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    x = "Observed attrition risk",
    y = "Predicted attrition risk"
  )

ggsave(
  "outputs/design/internal_validation/out/cv_obs_pred_risk_plot.png",
  cv_risk_plot,
  width = 6,
  height = 6,
  units = "in"
)

# Observed versus predicted counts ---------------------------------------------

cv_count_plot <- cv_df %>%
  ggplot(aes(x = r_obs, y = count_pred_mean)) +
  geom_abline(linetype = "dashed") +
  geom_errorbar(
    aes(ymin = count_pp_q025, ymax = count_pp_q975),
    alpha = 0.25,
    width = 0
  ) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(
    x = "Observed attrition count",
    y = "Predicted attrition count"
  )

ggsave(
  "outputs/design/internal_validation/out/cv_obs_pred_count_plot.png",
  cv_count_plot,
  width = 6,
  height = 6,
  units = "in"
)

# Count recovery by trial -------------------------------------------------------

cv_count_recovery <- cv_df %>%
  arrange(r_obs) %>%
  mutate(plot_row = row_number()) %>%
  ggplot(aes(x = plot_row, y = r_obs)) +
  geom_point(alpha = 0.7) +
  geom_point(aes(y = count_pred_mean), colour = "skyblue") +
  geom_errorbar(
    aes(ymin = count_pp_q025, ymax = count_pp_q975),
    colour = "skyblue",
    alpha = 0.6,
    width = 0.1
  ) +
  theme_classic() +
  labs(
    x = "Trial, ordered by observed attrition count",
    y = "Attrition count"
  )

ggsave(
  "outputs/design/internal_validation/out/cv_count_recovery.png",
  cv_count_recovery,
  width = 10,
  height = 4,
  units = "in"
)

# Residual checks ---------------------------------------------------------------

cv_resid_pred <- cv_df %>%
  mutate(risk_error = risk_pred_mean - risk_obs) %>%
  ggplot(aes(x = risk_pred_mean, y = risk_error)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_classic() +
  labs(
    x = "Predicted attrition risk",
    y = "Prediction error on risk scale"
  )

ggsave(
  "outputs/design/internal_validation/out/cv_residuals_by_pred_risk.png",
  cv_resid_pred,
  width = 6,
  height = 5,
  units = "in"
)

cv_resid_n <- cv_df %>%
  mutate(risk_error = risk_pred_mean - risk_obs) %>%
  ggplot(aes(x = N, y = risk_error)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_classic() +
  labs(
    x = "Trial arm size",
    y = "Prediction error on risk scale",
  )

ggsave(
  "outputs/design/internal_validation/out/cv_residuals_by_sample_size.png",
  cv_resid_n,
  width = 6,
  height = 5,
  units = "in"
)

# Count concordance and calibration summaries ----------------------------------

z_cv_draws <- "outputs/design/internal_validation/out/cv_prediction_draws.rds"
z_out_dir <- "outputs/design/internal_validation/out"

dir.create(z_out_dir, recursive = TRUE, showWarnings = FALSE)

# Use cv_draws from the current run when available. This fallback allows this
# section to be rerun interactively after cv_prediction_draws.rds has been saved.
if (!exists("cv_draws")) {
  if (!file.exists(z_cv_draws)) {
    stop(
      "Could not find held-out posterior prediction draws: ",
      z_cv_draws,
      "\nRun scripts/06b_design_validate_output.R first."
    )
  }
  cv_draws <- readRDS(z_cv_draws)
}
# Helper functions -------------------------------------------------------------

summ_vec <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(
      tibble(
        mean = NA_real_,
        median = NA_real_,
        sd = NA_real_,
        q025 = NA_real_,
        q975 = NA_real_,
        n_posterior_draws_used = 0L
      )
    )
  }
  
  data.frame(
    mean = mean(x),
    median = median(x),
    sd = sd(x),
    q025 = q025(x),
    q975 = q975(x),
    n_posterior_draws_used = length(x)
  )
}

fit_calibration_draw <- function(pred, y_binom) {
  
  pred <- pmin(pmax(pred, 1e-6), 1 - 1e-6)
  z <- qlogis(pred)
  
  out <- c(
    citl_intercept = NA_real_,
    calibration_intercept = NA_real_,
    calibration_slope = NA_real_,
    failed_citl = 1,
    failed_calibration = 1
  )
  
  fit_citl <- tryCatch(
    suppressWarnings(
      glm.fit(
        x = matrix(1, nrow = length(z), ncol = 1),
        y = y_binom,
        family = binomial(),
        offset = z
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit_citl) &&
      isTRUE(fit_citl$converged) &&
      all(is.finite(fit_citl$coefficients))) {
    out["citl_intercept"] <- fit_citl$coefficients[1]
    out["failed_citl"] <- 0
  }
  
  fit_cal <- tryCatch(
    suppressWarnings(
      glm.fit(
        x = cbind(intercept = 1, z = z),
        y = y_binom,
        family = binomial()
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit_cal) &&
      isTRUE(fit_cal$converged) &&
      all(is.finite(fit_cal$coefficients))) {
    out["calibration_intercept"] <- fit_cal$coefficients[1]
    out["calibration_slope"] <- fit_cal$coefficients[2]
    out["failed_calibration"] <- 0
  }
  
  out
}

reshape_prediction_draws <- function(dat) {
  
  dat <- dat[order(dat$row, dat$draw), ]
  
  draw_ids <- sort(unique(dat$draw))
  n_draw <- length(draw_ids)
  
  if (nrow(dat) %% n_draw != 0) {
    stop("Prediction draw object cannot be reshaped into draw x trial format.")
  }
  
  n_trial <- nrow(dat) / n_draw
  
  draw_mat <- matrix(dat$draw, nrow = n_draw, ncol = n_trial)
  
  if (!all(draw_mat == draw_ids)) {
    stop("Draw identifiers are not aligned within trials.")
  }
  
  list(
    draw_ids = draw_ids,
    risk = matrix(dat$risk_pred, nrow = n_draw, ncol = n_trial),
    count = matrix(dat$count_pred, nrow = n_draw, ncol = n_trial),
    fold = matrix(dat$fold, nrow = n_draw, ncol = n_trial)[1, ],
    row = matrix(dat$row, nrow = n_draw, ncol = n_trial)[1, ],
    N = matrix(dat$N, nrow = n_draw, ncol = n_trial)[1, ],
    r_obs = matrix(dat$r_obs, nrow = n_draw, ncol = n_trial)[1, ]
  )
}

# Trial-level count concordance statistic
# This estimates the probability that, for two trials with different observed
# attrition counts, the trial with more attrition has the higher expected count
count_c_stat <- function(pred, observed) {

  keep <- is.finite(pred) & is.finite(observed)
  pred <- pred[keep]
  observed <- observed[keep]

  observed_values <- sort(unique(observed))
  pred_below <- numeric()
  numerator <- 0
  denominator <- 0

  for (value in observed_values) {
    pred_value <- pred[observed == value]

    if (length(pred_below) > 0) {
      pred_diff <- outer(pred_value, pred_below, "-")

      numerator <- numerator +
        sum(pred_diff > 0) +
        0.5 * sum(pred_diff == 0)

      denominator <- denominator + length(pred_diff)
    }

    pred_below <- c(pred_below, pred_value)
  }

  if (denominator == 0) {
    return(NA_real_)
  }

  numerator / denominator
}

# Trial-level bootstrap for uncertainty around count concordance
bootstrap_count_c <- function(
    pred,
    observed,
    n_boot = 2000L,
    seed = 123L
) {

  set.seed(seed)
  n_trials <- length(pred)

  map_dbl(seq_len(n_boot), function(b) {
    idx <- sample.int(n_trials, size = n_trials, replace = TRUE)

    count_c_stat(
      pred = pred[idx],
      observed = observed[idx]
    )
  })
}

summarise_c_bootstrap <- function(x) {

  x <- x[is.finite(x)]

  tibble(
    bootstrap_mean = mean(x),
    bootstrap_median = median(x),
    bootstrap_sd = sd(x),
    bootstrap_q025 = q025(x),
    bootstrap_q975 = q975(x),
    n_bootstrap_samples = length(x)
  )
}

# Estimate concordance in recovering trial-level attrition counts
make_count_c_summary <- function(
    prediction,
    pred,
    observed,
    n_boot = 2000L,
    seed = 123L
) {

  count_c_boot <- bootstrap_count_c(
    pred = pred,
    observed = observed,
    n_boot = n_boot,
    seed = seed
  )

  bind_cols(
    tibble(
      prediction = prediction,
      metric = "trial_count_concordance",
      estimate = count_c_stat(
        pred = pred,
        observed = observed
      )
    ),
    summarise_c_bootstrap(count_c_boot),
    tibble(
      n_held_out_trials = length(observed),
      n_distinct_observed_counts = length(unique(observed))
    )
  )
}

make_count_c_fold_summary <- function(
    prediction,
    pred,
    observed,
    fold,
    n_boot = 2000L,
    seed = 123L
) {

  map_dfr(sort(unique(fold)), function(k) {

    idx <- which(fold == k)
    count_c_boot <- bootstrap_count_c(
      pred = pred[idx],
      observed = observed[idx],
      n_boot = n_boot,
      seed = seed + k
    )

    bind_cols(
      tibble(
        prediction = prediction,
        fold = k,
        metric = "trial_count_concordance",
        estimate = count_c_stat(
          pred = pred[idx],
          observed = observed[idx]
        )
      ),
      summarise_c_bootstrap(count_c_boot),
      tibble(
        n_held_out_trials = length(idx),
        n_distinct_observed_counts = length(unique(observed[idx]))
      )
    )
  })
}

make_calibration_summary <- function(prediction, pred_obj) {
  
  y_binom <- cbind(pred_obj$r_obs, pred_obj$N - pred_obj$r_obs)
  
  cal_draws <- t(vapply(
    seq_len(nrow(pred_obj$risk)),
    function(s) fit_calibration_draw(pred_obj$risk[s, ], y_binom),
    numeric(5)
  ))
  
  citl_failed <- sum(cal_draws[, "failed_citl"] == 1, na.rm = TRUE)
  cal_failed <- sum(cal_draws[, "failed_calibration"] == 1, na.rm = TRUE)
  
  rbind(
    cbind(
      data.frame(
        prediction = prediction,
        metric = "calibration_in_the_large_intercept"
      ),
      summ_vec(cal_draws[, "citl_intercept"]),
      data.frame(n_failed_model_fits = citl_failed)
    ),
    cbind(
      data.frame(
        prediction = prediction,
        metric = "calibration_intercept"
      ),
      summ_vec(cal_draws[, "calibration_intercept"]),
      data.frame(n_failed_model_fits = cal_failed)
    ),
    cbind(
      data.frame(
        prediction = prediction,
        metric = "calibration_slope"
      ),
      summ_vec(cal_draws[, "calibration_slope"]),
      data.frame(n_failed_model_fits = cal_failed)
    )
  )
}

# Generate supplementary summaries --------------------------------------------

# Posterior draws are retained for calibration summaries. Count concordance is
# calculated from posterior mean predictions, integrating over unknown trial
# effects for the marginal predictions.
prediction_types <- unique(cv_draws$prediction)
n_count_c_boot <- 2000L

cv_by_prediction <- map(
  prediction_types,
  function(prediction) {
    reshape_prediction_draws(cv_draws[cv_draws$prediction == prediction, ])
  }
) %>%
  set_names(prediction_types)

# Confirm predictions are aligned before comparing marginal and typical results
prediction_ref <- cv_by_prediction[[1]]

prediction_aligned <- map_lgl(
  cv_by_prediction,
  function(x) {
    identical(x$row, prediction_ref$row) &&
      identical(x$fold, prediction_ref$fold) &&
      identical(x$N, prediction_ref$N) &&
      identical(x$r_obs, prediction_ref$r_obs)
  }
)

if (!all(prediction_aligned)) {
  stop("Marginal and typical held-out predictions are not aligned.")
}

# Posterior mean expected attrition count for each held-out trial
count_predictions <- map(
  cv_by_prediction,
  function(x) colMeans(x$count, na.rm = TRUE)
)

# Fold-specific reference prediction using trial size and the aggregate
# attrition risk from the corresponding training set
count_predictions$size_only <- vapply(
  seq_along(prediction_ref$row),
  function(i) {
    train_idx <- prediction_ref$fold != prediction_ref$fold[i]
    train_risk <- sum(prediction_ref$r_obs[train_idx]) /
      sum(prediction_ref$N[train_idx])

    prediction_ref$N[i] * train_risk
  },
  numeric(1)
)

count_c_summary <- imap_dfr(
  count_predictions,
  ~ make_count_c_summary(
    prediction = .y,
    pred = .x,
    observed = prediction_ref$r_obs,
    n_boot = n_count_c_boot,
    seed = 2000L + which(names(count_predictions) == .y)
  )
)

count_c_fold_summary <- imap_dfr(
  count_predictions,
  ~ make_count_c_fold_summary(
    prediction = .y,
    pred = .x,
    observed = prediction_ref$r_obs,
    fold = prediction_ref$fold,
    n_boot = n_count_c_boot,
    seed = 3000L + 100L * which(names(count_predictions) == .y)
  )
)

calibration_summary <- map2_dfr(
  names(cv_by_prediction),
  cv_by_prediction,
  make_calibration_summary
)

count_c_summary <- count_c_summary[
  order(match(count_c_summary$prediction, names(count_predictions))),
]
count_c_fold_summary <- count_c_fold_summary[
  order(
    match(count_c_fold_summary$prediction, names(count_predictions)),
    count_c_fold_summary$fold
  ),
]
calibration_summary <- calibration_summary[
  order(calibration_summary$prediction, calibration_summary$metric),
]

# Save outputs -----------------------------------------------------------------

# Count concordance summaries
write_csv(
  count_c_summary,
  file.path(z_out_dir, "design_stage_cv_count_concordance.csv"),
  na = ""
)

saveRDS(
  count_c_summary,
  file.path(z_out_dir, "design_stage_cv_count_concordance.rds")
)

write_csv(
  count_c_fold_summary,
  file.path(z_out_dir, "design_stage_cv_count_concordance_by_fold.csv"),
  na = ""
)

saveRDS(
  count_c_fold_summary,
  file.path(z_out_dir, "design_stage_cv_count_concordance_by_fold.rds")
)

write_csv(
  calibration_summary,
  file.path(z_out_dir, "design_stage_cv_calibration_intercept_slope.csv"),
  na = ""
)

saveRDS(
  calibration_summary,
  file.path(z_out_dir, "design_stage_cv_calibration_intercept_slope.rds")
)

# Log --------------------------------------------------------------------------

log_lines <- c(
  "Design-stage cross-validation supplementary validation outputs",
  "",
  paste0("Input posterior draw file: ", z_cv_draws),
  paste0(
    "Primary count prediction: marginal posterior mean expected count ",
    "integrated over parameter uncertainty and unknown trial effects"
  ),
  paste0(
    "Count sensitivity prediction: typical posterior mean expected count ",
    "with unknown trial effects set to their central values"
  ),
  paste0(
    "Count ranking metric: trial-level count concordance with ",
    "a fold-specific size-only reference"
  ),
  paste0(
    "Count concordance uncertainty: trial-level bootstrap with ",
    n_count_c_boot,
    " samples"
  ),
  paste0("Posterior-draw predictions used for calibration: yes"),
  paste0("Prediction summaries included: ", paste(prediction_types, collapse = ", ")),
  paste0(
    "Held-out trials: ",
    length(unique(cv_by_prediction[[1]]$row))
  ),
  paste0(
    "Observed attrition count: ",
    sum(cv_by_prediction[[1]]$r_obs)
  ),
  paste0(
    "Total randomised: ",
    sum(cv_by_prediction[[1]]$N)
  ),
  paste0(
    "Calibration model failed fits: ",
    sum(calibration_summary$n_failed_model_fits, na.rm = TRUE)
  ),
  "",
  "Outputs produced:",
  file.path(z_out_dir, "design_stage_cv_count_concordance.csv"),
  file.path(z_out_dir, "design_stage_cv_count_concordance.rds"),
  file.path(z_out_dir, "design_stage_cv_count_concordance_by_fold.csv"),
  file.path(z_out_dir, "design_stage_cv_count_concordance_by_fold.rds"),
  file.path(z_out_dir, "design_stage_cv_calibration_intercept_slope.csv"),
  file.path(z_out_dir, "design_stage_cv_calibration_intercept_slope.rds"),
  "",
  paste(
    "Interpretation note: count concordance is the probability that the",
    "trial with the higher observed attrition count also has the higher",
    "posterior mean expected count, excluding pairs tied on observed count."
  )
)

writeLines(
  log_lines,
  file.path(z_out_dir, "design_stage_cv_supplementary_README.txt")
)
