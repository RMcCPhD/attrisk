
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

# Grouped-data C-statistic.
# This is an arm-level, participant-weighted risk-ranking metric:
# it estimates the probability that an observed attrition event is assigned
# a higher predicted trial-arm risk than an observed non-event.
grouped_c_stat <- function(pred, events, non_events) {
  
  keep <- is.finite(pred) &
    is.finite(events) &
    is.finite(non_events) &
    events >= 0 &
    non_events >= 0
  
  pred <- pred[keep]
  events <- events[keep]
  non_events <- non_events[keep]
  
  n_events <- sum(events)
  n_non_events <- sum(non_events)
  
  if (n_events <= 0 || n_non_events <= 0) {
    return(NA_real_)
  }
  
  ord <- order(pred)
  pred <- pred[ord]
  events <- events[ord]
  non_events <- non_events[ord]
  
  ties <- rle(pred)
  tie_end <- cumsum(ties$lengths)
  tie_start <- tie_end - ties$lengths + 1L
  
  numerator <- 0
  non_events_below <- 0
  
  for (g in seq_along(tie_start)) {
    idx <- tie_start[g]:tie_end[g]
    events_tie <- sum(events[idx])
    non_events_tie <- sum(non_events[idx])
    
    numerator <- numerator +
      events_tie * (non_events_below + 0.5 * non_events_tie)
    
    non_events_below <- non_events_below + non_events_tie
  }
  
  numerator / (n_events * n_non_events)
}

make_grouped_roc <- function(pred, events, non_events) {
  
  tibble(
    pred = pred,
    events = events,
    non_events = non_events
  ) %>%
    filter(
      is.finite(pred),
      is.finite(events),
      is.finite(non_events)
    ) %>%
    arrange(desc(pred)) %>%
    group_by(pred) %>%
    reframe(
      events = sum(events),
      non_events = sum(non_events)
    ) %>%
    mutate(
      sensitivity = cumsum(events) / sum(events),
      one_minus_specificity = cumsum(non_events) / sum(non_events)
    ) %>%
    bind_rows(
      tibble(
        pred = Inf,
        events = 0,
        non_events = 0,
        sensitivity = 0,
        one_minus_specificity = 0
      ),
      .,
      tibble(
        pred = -Inf,
        events = 0,
        non_events = 0,
        sensitivity = 1,
        one_minus_specificity = 1
      )
    )
}

# Draw-level posterior predictive performance
# Rows = posterior draws
# Columns = trial arms
r_pred_draws <- ppc_out$r_new

performance_draws <- map_dfr(seq_len(nrow(r_pred_draws)), function(s) {
  
  r_pred_s <- r_pred_draws[s, ]
  
  error <- r_pred_s - imp_dat$r
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
    cor_obs_pred = safe_cor(imp_dat$r, r_pred_s)
  )
})

# Arm-level discrimination from latent attrition probabilities -----------------

auc_draws <- map_dbl(seq_len(nrow(ppc_out$p)), function(s) {
  
  grouped_c_stat(
    pred = ppc_out$p[s, ],
    events = imp_dat$r,
    non_events = imp_dat$N - imp_dat$r
  )
  
})

auc_summary <- posterior_scalar_summary(auc_draws) %>%
  mutate(
    metric = "grouped_c_statistic_auc",
    n_posterior_draws_used = sum(is.finite(auc_draws)),
    n_arms = imp_dat$I,
    n_observed_attrition_events = sum(imp_dat$r),
    n_observed_non_events = sum(imp_dat$N - imp_dat$r),
    .before = 1
  )

write_csv(
  auc_summary,
  "outputs/design/posterior_auc.csv"
)

saveRDS(
  auc_summary,
  "outputs/design/posterior_auc.rds"
)

auc_curve <- make_grouped_roc(
  pred = ppc_tbl$p_mean,
  events = imp_dat$r,
  non_events = imp_dat$N - imp_dat$r
)

write_csv(
  auc_curve,
  "outputs/design/posterior_auc_curve.csv"
)

auc_plot <- auc_curve %>%
  ggplot(aes(x = one_minus_specificity, y = sensitivity)) +
  geom_abline(linetype = "dashed") +
  geom_step(linewidth = 0.7, colour = "skyblue") +
  coord_equal() +
  theme_classic() +
  annotate(
    "text",
    x = 0.65,
    y = 0.15,
    label = paste0("C-statistic = ", round(auc_summary$mean, 3))
  ) +
  labs(
    x = "1 - specificity",
    y = "Sensitivity"
  )

ggsave(
  "outputs/design/posterior_auc_plot.png",
  auc_plot,
  width = 5,
  height = 5,
  units = "in"
)

# Point estimate coverage from posterior predictive intervals
# This is not draw-specific because it is defined using the full posterior interval
coverage_point <- ppc_tbl %>%
  summarise(
    covered_95 = mean(r_obs >= r_pred_q025 & r_obs <= r_pred_q975)
  ) %>%
  pull(covered_95)

# Summarise posterior distribution of performance metrics
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
    auc_summary %>%
      select(metric, mean, median, sd, q025, q975)
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
