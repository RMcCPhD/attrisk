
library(tidyverse)
library(broom)
library(rstan)
library(bayesplot)
library(tidybayes)
library(posterior)
library(markdown)
library(flexsurv)
library(metafor)
library(multinma)

options(scipen = 999)
Sys.setenv(STANC3 = TRUE)

# Function to check dataset missingness quickly
check_nas <- function(dat) {
  dat %>% 
    reframe(across(everything(), ~ sum(is.na(.)))) %>% 
    pivot_longer(everything())
}

# Function to find the nearest neighbour IPD trial to an aggregate trial
# Use the fu_ratio of the nearest neighbour as an approximate fu ratio for the
# respective aggregate trials
nn_ratio <- function(row, ipd) {
  
  # Calculated squared Euclidean distances
  distances <- (ipd$scale_arm_n - row$scale_arm_n) ^ 2 +
               (ipd$scale_arm_attr - row$scale_arm_attr) ^ 2 +
               (ipd$scale_dur - row$scale_dur) ^ 2
  
  # Get match (smallest distance)
  best_match <- which.min(distances)
  
  # Return mft ratio from best match
  return(ipd$mft_ratio[best_match])
  
}

# Function to produce a tidy count of each variable, informing coding
# Contingency for dummy variables, summary for numeric
tidy_sum <- function(df, var_name) {
  var <- df[[var_name]]
  cat("Variable:", var_name, "\n")
  if (is.integer(var)) {
    tbl <- table(var) %>% as.data.frame()
    print(tbl)
  } else if (is.numeric(var)) {
    tbl <- summary(var)
    print(tbl)
  } else {NA}
}

# Function for covariate posterior predictive checks
plot_subgroup_ppcs <- function(y, yrep, covariate_df, covariate_names, stat = "mean") {
  stopifnot(length(y) == ncol(yrep))
  
  plots <- list()
  
  for (var in covariate_names) {
    group <- covariate_df[[var]]
    
    # Skip if group has too many levels
    if (length(unique(group)) > 10) {
      warning(glue::glue("Skipping {var} (too many levels)"))
      next
    }
    
    p <- ppc_stat_grouped(
      y = y,
      yrep = yrep,
      group = group,
      stat = stat
    ) +
      ggtitle(glue::glue("PPC by {var}")) +
      theme_minimal()
    
    plots[[var]] <- p
  }
  
  return(plots)
}

# Helper functions for internal cross-validation -------------------------------

# Gompertz CDF
# Numerical stability baked in, reduces to exponential if shape is near-zero
gompertz_prob <- function(b, a, dt) {
  
  ifelse(
    abs(a) < 1e-4,
    1 - exp(-b * dt),
    1 - exp(-b / a * (exp(a * dt) - 1))
  )
  
}

# Estimate lower quantile
q025 <- function(x) {
  
  quantile(x, 0.025, na.rm = TRUE, names = FALSE)
  
}

# Estimate upper quantile
q975 <- function(x) {
  
  quantile(x, 0.975, na.rm = TRUE, names = FALSE)
  
}

# Estimate lower quantile by column
col_q025 <- function(x) {
  
  apply(x, 2, q025)
  
}

# Estimate upper quantile by column
col_q975 <- function(x) {
  
  apply(x, 2, q975)
  
}

# Safely estimate correlation
safe_cor <- function(x, y) {
  
  ok <- complete.cases(x, y)
  
  if (sum(ok) < 2 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    
    return(NA_real_)
    
  }
  
  cor(x[ok], y[ok], method = "spearman")
  
}

# Reshape fold predictions to draw-level long format
make_draw_df <- function(
    prediction,
    fold,
    test_idx,
    N_test,
    r_test,
    count_exp_mat
) {
  
  S <- nrow(count_exp_mat)
  n_test <- ncol(count_exp_mat)
  
  tibble(
    prediction = prediction,
    fold = fold,
    draw = rep(seq_len(S), times = n_test),
    row = rep(test_idx, each = S),
    N = rep(N_test, each = S),
    r_obs = rep(r_test, each = S),
    count_pred = as.vector(count_exp_mat)
  ) %>%
    mutate(
      risk_obs = r_obs / N,
      risk_pred = count_pred / N
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
  fit_path <- file.path(cv_dir, sprintf("folds_models/fit_fold_%d.rds", k))
  
  if (!file.exists(fit_path)) {
    
    stop("Could not find fitted fold model: ", fit_path)
    
  }
  
  fit_k <- readRDS(fit_path)
  
  # Extract prediction parameters
  ext <- extract(fit_k, pars = c("mu_mu1", "beta"))
  S_all <- length(ext$mu_mu1)
  
  # Use all posterior draws unless max_draws is specified
  if (is.finite(max_draws) && max_draws < S_all) {
    
    draw_idx <- sample(seq_len(S_all), size = max_draws)
    
  } else {
    
    draw_idx <- seq_len(S_all)
    
  }
  
  # Define parameters
  mu_mu1 <- ext$mu_mu1[draw_idx]
  beta_arr <- ext$beta[draw_idx, , drop = FALSE]
  
  S <- length(draw_idx)
  n_test <- length(test_idx)
  
  # Held-out test data
  X_test <- imp_dat$X[test_idx, , drop = FALSE]
  dt_test <- imp_dat$dt[test_idx]
  N_test <- imp_dat$N[test_idx]
  r_test <- imp_dat$r[test_idx]
  
  # Shape-model parameters
  alpha_shape <- imp_dat$alpha_shape
  beta_shape_lograte <- imp_dat$beta_lograte
  
  # Store typical expected counts
  pred_count_exp_typical <- matrix(NA_real_, nrow = S, ncol = n_test)
  
  # Posterior typical-prediction loop
  for (s in seq_len(S)) {
    
    # Define typical log-rate and shape parameters
    mu1_typical <- rep(mu_mu1[s], n_test)
    
    log_b_typical <- as.vector(mu1_typical + X_test %*% beta_arr[s, ])
    b_typical <- exp(log_b_typical)
    
    a_typical <- alpha_shape + beta_shape_lograte * log_b_typical
    
    # Estimate attrition probability
    p_typical <- gompertz_prob(
      b = b_typical,
      a = a_typical,
      dt = dt_test
    )
    
    # Estimate typical expected counts
    pred_count_exp_typical[s, ] <- N_test * p_typical
    
  }
  
  # Row-level held-out prediction summaries
  pred_sum <- tibble(
    fold = k,
    row = test_idx,
    N = N_test,
    r_obs = r_test,
    dt = dt_test,
    
    # Typical expected counts
    count_typical_mean = colMeans(pred_count_exp_typical),
    count_typical_median = apply(pred_count_exp_typical, 2, median)
    
  ) %>%
    mutate(
      risk_obs = r_obs / N,
      risk_typical_mean = count_typical_mean / N,
      risk_typical_median = count_typical_median / N
    )
  
  # Draw-level predictions for credible intervals around summary metrics
  pred_draws <- make_draw_df(
    prediction = "typical",
    fold = k,
    test_idx = test_idx,
    N_test = N_test,
    r_test = r_test,
    count_exp_mat = pred_count_exp_typical
  )
  
  list(
    pred = pred_sum,
    draws = pred_draws
  )
  
}

# Summarise finite posterior draws
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

# Fit count calibration models within a posterior draw
fit_count_calibration_draw <- function(pred, observed) {
  
  pred <- pmax(pred, 1e-6)
  z <- log(pred)
  
  out <- c(
    count_citl_intercept = NA_real_,
    count_calibration_intercept = NA_real_,
    count_calibration_slope = NA_real_,
    failed_citl = 1,
    failed_calibration = 1
  )
  
  # Fit calibration-in-the-large model
  fit_citl <- tryCatch(
    suppressWarnings(
      glm.fit(
        x = matrix(1, nrow = length(z), ncol = 1),
        y = observed,
        family = poisson(),
        offset = z
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit_citl) &&
      isTRUE(fit_citl$converged) &&
      all(is.finite(fit_citl$coefficients))) {
    
    out["count_citl_intercept"] <- fit_citl$coefficients[1]
    out["failed_citl"] <- 0
    
  }
  
  # Fit calibration intercept and slope model
  fit_cal <- tryCatch(
    suppressWarnings(
      glm.fit(
        x = cbind(intercept = 1, z = z),
        y = observed,
        family = poisson()
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit_cal) &&
      isTRUE(fit_cal$converged) &&
      all(is.finite(fit_cal$coefficients))) {
    
    out["count_calibration_intercept"] <- fit_cal$coefficients[1]
    out["count_calibration_slope"] <- fit_cal$coefficients[2]
    out["failed_calibration"] <- 0
    
  }
  
  out
  
}

# Fit risk calibration models using binomial arm-level outcomes
# CITL: logit(observed risk) = intercept + offset(logit(predicted risk))
# Joint calibration: logit(observed risk) = intercept + slope * logit(predicted risk)
# The two-column response retains each arm's event count and denominator. Unlike
# equal-arm error summaries, these likelihood fits account for arm sample size
# Calibration coefficients are on the logit scale (ideal intercept 0, slope 1)
fit_risk_calibration_draw <- function(pred, observed, N) {
  
  out <- c(
    risk_citl_intercept = NA_real_,
    risk_calibration_intercept = NA_real_,
    risk_calibration_slope = NA_real_,
    failed_citl = 1,
    failed_calibration = 1
  )
  
  if (length(pred) != length(observed) || length(N) != length(observed) ||
      length(pred) == 0 || any(!is.finite(pred)) ||
      any(!is.finite(observed)) || any(!is.finite(N)) ||
      any(N <= 0 | observed < 0 | observed > N) ||
      any(pred < 0 | pred > 1) ||
      sum(observed) == 0 || sum(N - observed) == 0) {
    
    return(out)
    
  }
  
  # Clip only for taking logits; saved risks and error metrics are not clipped
  z <- qlogis(pmin(pmax(pred, 1e-6), 1 - 1e-6))
  y <- cbind(events = observed, non_events = N - observed)
  
  fit_citl <- tryCatch(
    suppressWarnings(
      glm.fit(
        x = matrix(1, nrow = length(z), ncol = 1),
        y = y,
        family = binomial(link = "logit"),
        offset = z,
        control = list(maxit = 100)
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit_citl) && isTRUE(fit_citl$converged) &&
      !isTRUE(fit_citl$boundary) &&
      all(is.finite(fit_citl$coefficients))) {
    
    out["risk_citl_intercept"] <- fit_citl$coefficients[1]
    out["failed_citl"] <- 0
    
  }
  
  fit_cal <- tryCatch(
    suppressWarnings(
      glm.fit(
        x = cbind(intercept = 1, z = z),
        y = y,
        family = binomial(link = "logit"),
        control = list(maxit = 100)
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit_cal) && isTRUE(fit_cal$converged) &&
      !isTRUE(fit_cal$boundary) && fit_cal$rank == 2 &&
      all(is.finite(fit_cal$coefficients))) {
    
    out["risk_calibration_intercept"] <- fit_cal$coefficients[1]
    out["risk_calibration_slope"] <- fit_cal$coefficients[2]
    out["failed_calibration"] <- 0
    
  }
  
  out
  
}

# Reshape posterior draws to draw-by-trial matrices
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
    count = matrix(dat$count_pred, nrow = n_draw, ncol = n_trial),
    risk = matrix(dat$count_pred / dat$N, nrow = n_draw, ncol = n_trial),
    fold = matrix(dat$fold, nrow = n_draw, ncol = n_trial)[1, ],
    row = matrix(dat$row, nrow = n_draw, ncol = n_trial)[1, ],
    N = matrix(dat$N, nrow = n_draw, ncol = n_trial)[1, ],
    r_obs = matrix(dat$r_obs, nrow = n_draw, ncol = n_trial)[1, ],
    risk_obs = matrix(dat$r_obs / dat$N, nrow = n_draw, ncol = n_trial)[1, ]
  )
  
}

# Trial-level concordance on either scale
# Observed ties are excluded; predicted ties receive half credit. Each comparable
# pair of arms has equal weight. Risk concordance ranks observed r/N against
# predicted risk
# Comparable pairs are prepared once per analysis to reuse across posterior draws
make_concordance_pairs <- function(observed) {
  
  which(outer(observed, observed, ">"), arr.ind = TRUE)
  
}

trial_c_stat <- function(pred, observed, pairs = NULL, pair_weights = NULL) {
  
  if (length(pred) != length(observed)) {
    
    stop("Predicted and observed values must have the same length.")
    
  }
  
  if (is.null(pairs)) {
    
    pairs <- make_concordance_pairs(observed)
    
  }
  
  if (nrow(pairs) == 0) {
    
    return(NA_real_)
    
  }
  
  difference <- pred[pairs[, 1]] - pred[pairs[, 2]]
  if (is.null(pair_weights)) pair_weights <- rep(1, nrow(pairs))
  stopifnot(
    length(pair_weights) == nrow(pairs),
    all(is.finite(pair_weights)),
    all(pair_weights > 0)
  )
  valid <- is.finite(difference)
  difference <- difference[valid]
  pair_weights <- pair_weights[valid]
  if (length(difference) == 0) {
    
    return(NA_real_)
    
  }
  
  weighted.mean((difference > 0) + 0.5 * (difference == 0), pair_weights)
  
}

# Summarise the posterior distribution of concordance
summarise_c_posterior <- function(x) {
  
  s <- summ_vec(x)
  tibble(
    posterior_mean = s$mean,
    posterior_median = s$median,
    posterior_sd = s$sd,
    posterior_q025 = s$q025,
    posterior_q975 = s$q975,
    n_posterior_draws = s$n_posterior_draws_used
  )
  
}

# Retain the original count metadata name and use its risk counterpart as needed
label_distinct_outcome <- function(dat, scale) {
  
  names(dat)[names(dat) == "n_distinct_observed"] <-
    paste0(
      "n_distinct_observed_",
      if (scale == "count") "counts" else "risks"
    )
  dat
  
}

# Estimate posterior concordance, selecting predictions and observations together
make_c_summary <- function(prediction, pred_obj, scale = c("count", "risk"),
                           size_weighted = FALSE) {
  
  scale <- match.arg(scale)
  pred_mat <- pred_obj[[scale]]
  observed <- if (scale == "count") pred_obj$r_obs else pred_obj$risk_obs
  pairs <- make_concordance_pairs(observed)
  pair_weights <- NULL
  if (size_weighted) {
    
    stopifnot(
      length(pred_obj$N) == length(observed),
      all(is.finite(pred_obj$N)),
      all(pred_obj$N > 0)
    )
    pair_weights <- as.numeric(pred_obj$N[pairs[, 1]]) *
      as.numeric(pred_obj$N[pairs[, 2]])
    
  }
  
  c_draws <- map_dbl(
    seq_len(nrow(pred_mat)),
    function(s) {
      
      trial_c_stat(pred_mat[s, ], observed, pairs, pair_weights)
      
    }
  )
  
  bind_cols(
    tibble(
      prediction = prediction,
      metric = paste0(
        "trial_",
        scale,
        "_concordance",
        if (size_weighted) "_size_weighted" else ""
      ),
      estimate = trial_c_stat(
        colMeans(pred_mat),
        observed,
        pairs,
        pair_weights
      )
    ),
    summarise_c_posterior(c_draws),
    tibble(
      n_held_out_trials = length(observed),
      n_distinct_observed = length(unique(observed))
    )
  ) %>%
    label_distinct_outcome(scale)
  
}

# Subset the arm dimension without dropping single-arm matrices
subset_prediction_object <- function(pred_obj, idx) {
  
  list(
    count = pred_obj$count[, idx, drop = FALSE],
    risk = pred_obj$risk[, idx, drop = FALSE],
    N = pred_obj$N[idx],
    r_obs = pred_obj$r_obs[idx],
    risk_obs = pred_obj$risk_obs[idx],
    fold = pred_obj$fold[idx]
  )
  
}

make_c_fold_summary <- function(prediction, pred_obj, scale = c(
  "count",
  "risk"
),
size_weighted = FALSE) {
  
  scale <- match.arg(scale)
  map_dfr(
    sort(unique(pred_obj$fold)),
    function(k) {
      
      idx <- which(pred_obj$fold == k)
      make_c_summary(
        prediction,
        subset_prediction_object(pred_obj, idx),
        scale,
        size_weighted = size_weighted
      ) %>%
        mutate(fold = k, .after = prediction)
      
    }
  )
  
}

# Summarise calibration fits across posterior prediction draws
# Count calibration keeps the original Poisson/log specification; risk
# calibration uses binomial/logit models and the original arm denominators
make_calibration_summary <- function(
    prediction, pred_obj, scale = c("count", "risk")
) {
  
  scale <- match.arg(scale)
  pred_mat <- pred_obj[[scale]]
  cal_draws <- t(vapply(
    seq_len(nrow(pred_mat)),
    function(s) {
      
      if (scale == "count") {
        
        fit_count_calibration_draw(pred_mat[s, ], pred_obj$r_obs)
        
      } else {
        
        fit_risk_calibration_draw(pred_mat[s, ], pred_obj$r_obs, pred_obj$N)
        
      }
      
    },
    numeric(5)
  ))
  
  citl_failed <- sum(cal_draws[, "failed_citl"] == 1, na.rm = TRUE)
  cal_failed <- sum(cal_draws[, "failed_calibration"] == 1, na.rm = TRUE)
  
  bind_rows(
    bind_cols(
      tibble(
        prediction = prediction,
        metric = paste0(scale, "_calibration_in_the_large_intercept")
      ),
      summ_vec(cal_draws[, paste0(scale, "_citl_intercept")]),
      tibble(n_failed_model_fits = citl_failed)
    ),
    bind_cols(
      tibble(
        prediction = prediction,
        metric = paste0(scale, "_calibration_intercept")
      ),
      summ_vec(cal_draws[, paste0(scale, "_calibration_intercept")]),
      tibble(n_failed_model_fits = cal_failed)
    ),
    bind_cols(
      tibble(
        prediction = prediction,
        metric = paste0(scale, "_calibration_slope")
      ),
      summ_vec(cal_draws[, paste0(scale, "_calibration_slope")]),
      tibble(n_failed_model_fits = cal_failed)
    )
  )
  
}

make_calibration_fold_summary <- function(
    prediction, pred_obj, scale = c("count", "risk")
) {
  
  scale <- match.arg(scale)
  map_dfr(
    sort(unique(pred_obj$fold)),
    function(k) {
      
      idx <- which(pred_obj$fold == k)
      make_calibration_summary(
        prediction,
        subset_prediction_object(pred_obj, idx),
        scale
      ) %>%
        mutate(fold = k, n_held_out_trials = length(idx), .after = prediction)
      
    }
  )
  
}

# Distribution parameter helpers -----------------------------------------------

get_lnorm_meanlog <- function(mu, d) {log(mu^2 / sqrt(d^2 + mu^2))}
get_lnorm_sdlog <- function(mu, d) {sqrt(log(1 + d^2 / mu^2))}
get_gamma_shape <- function(mu, d) {(mu / d)^2}
get_gamma_rate <- function(mu, d) {mu / d^2}

get_weibull_shape <- function(mu, d) {
  cv <- d / mu
  map_dbl(
    cv, ~ {
      uniroot(
        function(k) {sqrt(gamma(1 + 2 / k) / gamma(1 + 1 / k)^2 - 1) - .x},
        interval = c(0.1, 50)
      )$root
    }
  )
}

get_weibull_scale <- function(mu, shape) {mu / gamma(1 + 1 / shape)}
