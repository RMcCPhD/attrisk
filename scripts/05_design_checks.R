
# Posterior predictive checks for full design stage model
# Trace plots not generated because Rhat = 1 achieved for all parameters
# Assessing:
# - Recovery of attrition counts
# - Apparent count performance
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

# Posterior expected attrition counts
count_pred_draws <- sweep(ppc_out$p, 2, imp_dat$N, "*")

# Helper for row-level posterior summaries
# For r_new these are posterior predictive summaries
# For p these are posterior summaries of the latent attrition probability
posterior_col_summary <- function(x, prefix) {
  tibble(
    "{prefix}_mean" := colMeans(x, na.rm = TRUE),
    "{prefix}_median" := apply(x, 2, median, na.rm = TRUE),
    "{prefix}_sd" := apply(x, 2, sd, na.rm = TRUE),
    "{prefix}_q025" := apply(
      x,
      2,
      quantile,
      probs = 0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    "{prefix}_q975" := apply(
      x,
      2,
      quantile,
      probs = 0.975,
      na.rm = TRUE,
      names = FALSE
    )
  )
}

ppc_tbl <- tibble(
  row = seq_len(imp_dat$I),
  r_obs = imp_dat$r,
  N = imp_dat$N,
  obs_risk = r_obs / N
) %>%
  bind_cols(
    posterior_col_summary(ppc_out$r_new, "r_pred"),
    posterior_col_summary(count_pred_draws, "count_pred"),
    posterior_col_summary(ppc_out$p, "p")
  )

write_csv(
  ppc_tbl %>% mutate(across(everything(), ~ round(., 3))), 
  "outputs/design/posterior_summary.csv"
)
# Performance summary
# Helper for posterior summaries of scalar performance metrics
posterior_scalar_summary <- function(x) {
  tibble(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    q975 = quantile(x, 0.975, na.rm = TRUE, names = FALSE)
  )
}

# Helper for safe correlations by posterior draw
safe_cor <- function(x, y) {
  ok <- complete.cases(x, y)
  
  if (sum(ok) < 2 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(NA_real_)
  }
  
  cor(x[ok], y[ok])
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

  tibble(
    prediction = prediction,
    metric = "trial_count_concordance",
    estimate = count_c_stat(
      pred = pred,
      observed = observed
    ),
    bootstrap_mean = mean(count_c_boot, na.rm = TRUE),
    bootstrap_median = median(count_c_boot, na.rm = TRUE),
    bootstrap_sd = sd(count_c_boot, na.rm = TRUE),
    bootstrap_q025 = quantile(
      count_c_boot,
      0.025,
      na.rm = TRUE,
      names = FALSE
    ),
    bootstrap_q975 = quantile(
      count_c_boot,
      0.975,
      na.rm = TRUE,
      names = FALSE
    ),
    n_bootstrap_samples = sum(is.finite(count_c_boot)),
    n_trials = length(observed),
    n_distinct_observed_counts = n_distinct(observed)
  )
}

# Draw-level posterior expected count performance
# Rows = posterior draws
# Columns = trials

performance_draws <- map_dfr(seq_len(nrow(count_pred_draws)), function(s) {
  
  count_pred_s <- count_pred_draws[s, ]
  
  error <- count_pred_s - imp_dat$r
  abs_error <- abs(error)
  sq_error <- error^2
  rel_error <- error / pmax(imp_dat$r, 1)
  abs_rel_error <- abs(rel_error)
  
  tibble(
    draw = s,
    n = imp_dat$I,
    mean_error = mean(error),
    rel_error = mean(rel_error),
    mae = mean(abs_error),
    rmse = sqrt(mean(sq_error)),
    cor_obs_pred = safe_cor(imp_dat$r, count_pred_s)
  )
})

# Trial-level count concordance ------------------------------------------------

# Fitted prediction includes trial effects informed by the observed count
# This is an apparent development-data check rather than new-trial performance
count_predictions <- list(
  fitted_apparent = colMeans(count_pred_draws, na.rm = TRUE),
  size_only = imp_dat$N * sum(imp_dat$r) / sum(imp_dat$N)
)

n_count_c_boot <- 2000L

count_c_summary <- imap_dfr(
  count_predictions,
  ~ make_count_c_summary(
    prediction = .y,
    pred = .x,
    observed = imp_dat$r,
    n_boot = n_count_c_boot,
    seed = 123L + which(names(count_predictions) == .y)
  )
)

write_csv(
  count_c_summary,
  "outputs/design/posterior_count_concordance.csv"
)

saveRDS(
  count_c_summary,
  "outputs/design/posterior_count_concordance.rds"
)

# Point estimate coverage from posterior predictive intervals
# This is not draw-specific because it is defined using the full posterior interval
coverage_point <- ppc_tbl %>%
  summarise(
    covered_95 = mean(r_obs >= r_pred_q025 & r_obs <= r_pred_q975)
  ) %>%
  pull(covered_95)

# Summarise posterior distribution of performance metrics
# Concordance uncertainty is from the trial-level bootstrap
ppc_performance <- performance_draws %>%
  select(-draw) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  reframe(posterior_scalar_summary(value)) %>%
  bind_rows(
    tibble(
      metric = "covered_95",
      mean = coverage_point,
      median = NA_real_,
      sd = NA_real_,
      q025 = NA_real_,
      q975 = NA_real_
    )
  ) %>%
  bind_rows(
    count_c_summary %>%
      filter(prediction == "fitted_apparent") %>%
      transmute(
        metric = "trial_count_concordance_apparent",
        mean = estimate,
        median = bootstrap_median,
        sd = bootstrap_sd,
        q025 = bootstrap_q025,
        q975 = bootstrap_q975
      )
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write_csv(
  ppc_performance,
  "outputs/design/posterior_performance.csv"
)
# Plot posterior recovery
# Confirms adequate recovery
ppc_plot <- ppc_tbl %>% 
  ggplot(aes(x = row, y = r_obs)) +
  geom_point(size = 0.5) +
  geom_point(
    aes(x = row, y = count_pred_mean),
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
    x = "Trial",
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
