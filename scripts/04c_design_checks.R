
# Posterior predictive checks for design-stage model
# - Recovery of attrition counts and risks
# - Apparent count- and risk-scale performance in the development data
# - Covariate effect estimates (not PPC but here for initial check)

# Risk = count / control-arm sample size
# Error summaries give each control arm equal weight
# Aggregate risk summaries weight arms by size
# Expected predictions are used for performance metrics
# Posterior used for recovery intervals
# These are apparent performance metrics because outcome informs fit

rm(list = ls())
gc()

# Packages
library(tidyverse)
library(rstan)
library(bayesplot)
library(posterior)

# Define output directory
z_out_dir <- "outputs/design_stage/model/development_performance"
dir.create(z_out_dir, recursive = TRUE, showWarnings = FALSE)

# Import Stan data (must run 04b_design_fit.R within HPC)
imp_dat <- readRDS("processed_data/stan_data.rds")

# Check denominators are valid
if (length(imp_dat$N) != imp_dat$I ||
    length(imp_dat$r) != imp_dat$I ||
    any(!is.finite(imp_dat$N)) || any(imp_dat$N <= 0) ||
    any(imp_dat$N != floor(imp_dat$N)) ||
    any(!is.finite(imp_dat$r)) ||
    any(imp_dat$r < 0 | imp_dat$r > imp_dat$N) ||
    any(imp_dat$r != floor(imp_dat$r))) {
  stop(
    "Stan data must have positive N and r in [0, N] for every arm"
  )
}

# Import fitted model (cannot run without fitted model)
imp_mdl <- readRDS("outputs/design_stage/model/design_model.rds")

# Get posterior summaries
ppc_out <- extract(
  imp_mdl,
  pars = c("r_new", "p", "log_lik", "beta")
)

# Posterior expected predictions
risk_pred_draws <- ppc_out$p
count_pred_draws <- sweep(risk_pred_draws, 2, imp_dat$N, "*")

# Posterior predictive replicates
count_pp_draws <- ppc_out$r_new
risk_pp_draws <- sweep(count_pp_draws, 2, imp_dat$N, "/")

# List posterior
prediction_matrices <- list(
  risk_pred_draws = risk_pred_draws,
  count_pred_draws = count_pred_draws,
  count_pp_draws = count_pp_draws,
  risk_pp_draws = risk_pp_draws
)

# Helper for arm-level posterior summaries
posterior_col_summary <- function(x, prefix) {
  tibble(
    "{prefix}_mean" := colMeans(x, na.rm = TRUE),
    "{prefix}_median" := apply(x, 2, median, na.rm = TRUE),
    "{prefix}_sd" := apply(x, 2, sd, na.rm = TRUE),
    "{prefix}_q025" := apply(
      x, 2, quantile, probs = 0.025, na.rm = TRUE, names = FALSE
    ),
    "{prefix}_q975" := apply(
      x, 2, quantile, probs = 0.975, na.rm = TRUE, names = FALSE
    )
  )
}

# Get arm-level predictions
ppc_tbl <- tibble(
  row = seq_len(imp_dat$I),
  r_obs = imp_dat$r,
  N = imp_dat$N,
  obs_risk = r_obs / N,
  risk_obs = obs_risk
)

ppc_tbl <- ppc_tbl %>% 
  bind_cols(
    posterior_col_summary(count_pp_draws, "r_pred"),
    posterior_col_summary(count_pred_draws, "count_pred"),
    posterior_col_summary(risk_pred_draws, "p"),
    posterior_col_summary(risk_pred_draws, "risk_pred"),
    posterior_col_summary(risk_pp_draws, "risk_pp")
  )

# Save
write_csv(
  ppc_tbl, 
  file.path(z_out_dir, "posterior_level_counts_risks_sum.csv")
)

# Helpers for performance metrics ----------------------------------------------

# Posterior summary of a scalar metric
posterior_scalar_summary <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(tibble(
      mean = NA_real_, median = NA_real_, sd = NA_real_,
      q025 = NA_real_, q975 = NA_real_
    ))
  }
  tibble(
    mean = mean(x),
    median = median(x),
    sd = sd(x),
    q025 = quantile(x, 0.025, names = FALSE),
    q975 = quantile(x, 0.975, names = FALSE)
  )
}

# Safe Spearman correlation by posterior draw
safe_cor <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 2 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(NA_real_)
  }
  cor(x[ok], y[ok], method = "spearman")
}

# Control arm count and risk concordance
# Observed ties excluded, predicted ties given 0.5
# Two versions: Equal weighting, sample size-weighted
make_concordance_pairs <- function(observed) {
  which(outer(observed, observed, ">"), arr.ind = TRUE)
}

trial_c_stat <- function(pred, observed, pairs = NULL, pair_weights = NULL) {
  
  # Safeguards
  if (length(pred) != length(observed)) {
    stop("Predicted and observed must be the same length")
  }
  
  if (is.null(pairs)) {
    pairs <- make_concordance_pairs(observed)
  }
  
  if (nrow(pairs) == 0) {
    return(NA_real_)
  }
  
  # Calculate pair difference
  difference <- pred[pairs[, 1]] - pred[pairs[, 2]]
  if (is.null(pair_weights)) pair_weights <- rep(1, nrow(pairs))
  stopifnot(
    length(pair_weights) == nrow(pairs),
    all(is.finite(pair_weights)), all(pair_weights > 0)
  )
  
  valid <- is.finite(difference)
  difference <- difference[valid]
  pair_weights <- pair_weights[valid]
  
  if (length(difference) == 0) {
    return(NA_real_)
  }
  
  weighted.mean((difference > 0) + 0.5 * (difference == 0), pair_weights)
  
}

risk_obs <- imp_dat$r / imp_dat$N
count_pairs <- make_concordance_pairs(imp_dat$r)
risk_pairs <- make_concordance_pairs(risk_obs)

# Sample size weights for weighted concordance
risk_pair_weights <- as.numeric(imp_dat$N[risk_pairs[, 1]]) *
  as.numeric(imp_dat$N[risk_pairs[, 2]])

# Draw-level expected-prediction performance for counts and risks
# Constraints as safeguards - max(r_obs, 1); max(risk_obs, 1/n)
performance_draws <- map_dfr(seq_len(nrow(count_pred_draws)), function(s) {
  
  count_pred_s <- count_pred_draws[s, ]
  risk_pred_s <- risk_pred_draws[s, ]
  
  count_error <- count_pred_s - imp_dat$r
  risk_error <- risk_pred_s - risk_obs
  rel_count_error <- count_error / pmax(imp_dat$r, 1)
  rel_risk_error <- risk_error / pmax(risk_obs, 1 / imp_dat$N)
  
  tibble(
    draw = s,
    n = imp_dat$I,
    
    mean_obs_count = mean(imp_dat$r),
    mean_pred_count = mean(count_pred_s),
    total_obs_count = sum(imp_dat$r),
    total_pred_count = sum(count_pred_s),
    total_count_error = sum(count_error),
    mpe_count = mean(count_error),
    mae_count = mean(abs(count_error)),
    rmse_count = sqrt(mean(count_error^2)),
    mean_rel_error_count = mean(rel_count_error),
    mean_abs_rel_error_count = mean(abs(rel_count_error)),
    median_abs_rel_error_count = median(abs(rel_count_error)),
    cor_count = safe_cor(imp_dat$r, count_pred_s),
    trial_count_concordance_apparent = trial_c_stat(
      count_pred_s, imp_dat$r, count_pairs
    ),
    
    # Equal weight per control arm
    mean_obs_risk = mean(risk_obs),
    mean_pred_risk = mean(risk_pred_s),
    mpe_risk = mean(risk_error),
    mae_risk = mean(abs(risk_error)),
    rmse_risk = sqrt(mean(risk_error^2)),
    mean_rel_error_risk = mean(rel_risk_error),
    mean_abs_rel_error_risk = mean(abs(rel_risk_error)),
    median_abs_rel_error_risk = median(abs(rel_risk_error)),
    cor_risk = safe_cor(risk_obs, risk_pred_s),

    trial_risk_concordance_apparent = trial_c_stat(
      risk_pred_s, risk_obs, risk_pairs
    ),
    trial_risk_concordance_size_weighted_apparent = trial_c_stat(
      risk_pred_s, risk_obs, risk_pairs, risk_pair_weights
    ),
    
    # Sample-size-weighted
    aggregate_obs_risk = sum(imp_dat$r) / sum(imp_dat$N),
    aggregate_pred_risk = sum(count_pred_s) / sum(imp_dat$N),
    aggregate_risk_error = sum(count_error) / sum(imp_dat$N)
  )
  
})

# Get posterior draw performance
ppc_performance <- performance_draws %>%
  select(-draw) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  reframe(posterior_scalar_summary(value))

# Concordance summaries -------------------------------------------------------

# Create summary
make_c_summary <- function(prediction, pred_draws, observed, scale, 
                           pairs = NULL, pair_weights = NULL) {
  
  if (is.null(pairs)) {
    pairs <- make_concordance_pairs(observed)
  }
  
  c_draws <- map_dbl(seq_len(nrow(pred_draws)), function(s) {
    trial_c_stat(pred_draws[s, ], observed, pairs, pair_weights)
  })
  
  bind_cols(
    tibble(
      prediction = prediction,
      metric = paste0(
        "trial_", scale, "_concordance",
        if (!is.null(pair_weights)) "_size_weighted" else ""
      ),
      estimate = trial_c_stat(
        colMeans(pred_draws), 
        observed, pairs, pair_weights
      )
    ),
    posterior_scalar_summary(c_draws) %>% rename_with(~ paste0("posterior_", .x)),
    tibble(
      n_posterior_draws = sum(is.finite(c_draws)),
      n_trials = length(observed),
      n_distinct_observed = n_distinct(observed)
    )
  )
  
}

# Tibble to populate concordance summaries
make_c_reference_summary <- function(prediction, pred, observed, scale) {
  tibble(
    prediction = prediction,
    metric = paste0("trial_", scale, "_concordance"),
    estimate = trial_c_stat(pred, observed),
    posterior_mean = NA_real_,
    posterior_median = NA_real_,
    posterior_sd = NA_real_,
    posterior_q025 = NA_real_,
    posterior_q975 = NA_real_,
    n_posterior_draws = NA_integer_,
    n_trials = length(observed),
    n_distinct_observed = n_distinct(observed)
  )
}

aggregate_obs_risk <- sum(imp_dat$r) / sum(imp_dat$N)
size_only_prediction <- imp_dat$N * aggregate_obs_risk
aggregate_risk_prediction <- rep(aggregate_obs_risk, imp_dat$I)

count_c_summary <- bind_rows(
  make_c_summary(
    "fitted_apparent", count_pred_draws, imp_dat$r, "count", count_pairs
  ),
  make_c_reference_summary(
    "size_only", size_only_prediction, imp_dat$r, "count"
  )
) %>%
  rename(n_distinct_observed_counts = n_distinct_observed)

risk_c_summary <- bind_rows(
  make_c_summary(
    "fitted_apparent", risk_pred_draws, risk_obs, "risk", risk_pairs
  ),
  make_c_reference_summary(
    "aggregate_risk_only", aggregate_risk_prediction, risk_obs, "risk"
  )
) %>%
  rename(n_distinct_observed_risks = n_distinct_observed)

risk_c_size_weighted <- make_c_summary(
  "fitted_apparent", risk_pred_draws, risk_obs, 
  "risk", risk_pairs, risk_pair_weights
) %>%
  rename(n_distinct_observed_risks = n_distinct_observed)

# Posterior predictive interval coverage --------------------------------------

# Summarise coverage
interval_coverage <- bind_rows(
  ppc_tbl %>%
    transmute(
      scale = "count", observed = r_obs,
      lower = r_pred_q025, upper = r_pred_q975
    ),
  ppc_tbl %>%
    transmute(
      scale = "risk", observed = risk_obs,
      lower = risk_pp_q025, upper = risk_pp_q975
    )
) %>%
  group_by(scale) %>%
  summarise(
    n_trials = n(),
    n_covered_95 = sum(observed >= lower & observed <= upper),
    coverage_95 = mean(observed >= lower & observed <= upper),
    mean_interval_width = mean(upper - lower),
    .groups = "drop"
  )

# Combine performance summaries
ppc_performance <- bind_rows(
  ppc_performance,
  interval_coverage %>%
    transmute(
      metric = paste0("coverage_95_", scale),
      mean = coverage_95,
      median = NA_real_,
      sd = NA_real_,
      q025 = NA_real_,
      q975 = NA_real_
    )
)

write_csv(
  ppc_performance, 
  file.path(z_out_dir, "posterior_performance_sum.csv")
)

# Posterior risk recovery plot with coverage annotation
coverage_count <- filter(interval_coverage, scale == "count")
ppc_lbl <- sprintf(
  "PPI coverage: %d/%d control arms (%.1f%%)",
  coverage_count$n_covered_95,
  coverage_count$n_trials,
  100 * coverage_count$coverage_95
)

ppc_risk_plot <- ppc_tbl %>%
  arrange(risk_pred_mean) %>%
  mutate(plot_row = row_number()) %>%
  ggplot(aes(x = plot_row, y = risk_obs)) +
  geom_linerange(
    aes(ymin = risk_pp_q025, ymax = risk_pp_q975),
    colour = "skyblue", alpha = 0.7
  ) +
  geom_point(
    aes(y = risk_pred_mean), 
    colour = "skyblue", 
    size = 0.5
  ) +
  geom_point(size = 0.5) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1)
  ) +
  theme_classic() +
  labs(
    x = "Control arm, ordered by predicted risk", 
    y = "Attrition risk"
  ) +
  annotate(
    "text", x = Inf, y = Inf, label = ppc_lbl,
    hjust = 1.05, vjust = 1.5
  )

ggsave(
  file.path(z_out_dir, "posterior_risk_recovery.png"),
  ppc_risk_plot, 
  width = 10, 
  height = 4, 
  units = "in"
)
