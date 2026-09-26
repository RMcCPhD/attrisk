
# Generate typical outputs from internal cross-validation
# Performance, discrimination, and calibration summaries on count and risk scales
# Risk is evaluated at each control arm's follow-up (risk = count / sample size)

# Source configuration and functions
source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Set seed
set.seed(123)

# Define file directory for folds
cv_dir <- "hpc/design_cross_val/"
z_out_dir <- "outputs/design_stage/model/cross_validation/"

# Import Stan data
imp_stan <- readRDS("processed_data/stan_data.rds")

# Check arm counts and denominators
if (length(imp_stan$N) != imp_stan$I ||
    length(imp_stan$r) != imp_stan$I ||
    any(!is.finite(imp_stan$N)) || any(imp_stan$N <= 0) ||
    any(imp_stan$N != floor(imp_stan$N)) ||
    any(!is.finite(imp_stan$r)) ||
    any(imp_stan$r < 0 | imp_stan$r > imp_stan$N) ||
    any(imp_stan$r != floor(imp_stan$r))) {

  stop("Stan data must contain a positive integer N and an integer r in [0, N] for every arm.")

}

# Import fold metadata
imp_fold_ids <- readRDS(paste0(cv_dir, "folds/fold_id.rds"))

# Infer K from the saved allocation so all fitted folds are processed
fold_values <- sort(unique(imp_fold_ids))

# Check fold metadata
if (length(imp_fold_ids) != imp_stan$I) {

  stop("Fold IDs do not match the number of trials in the Stan data.")

}

if (anyNA(imp_fold_ids) ||
    !setequal(fold_values, seq_len(length(fold_values)))) {

  stop("Fold IDs must be contiguous integers starting at 1.")

}

K_folds <- length(fold_values)

# Import cv diagnostics
imp_cv_diag <- read_csv(paste0(cv_dir, "folds/cv_diagnostics.csv"))

# Generate cross-validation predictions ----------------------------------------

# Generate predictions for each fold
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

# Combine fold predictions
cv_df <- map_dfr(cv_out, "pred")
cv_draws <- map_dfr(cv_out, "draws")

# Save predictions
saveRDS(cv_df, paste0(z_out_dir, "cv_predictions.rds"))
saveRDS(cv_draws, paste0(z_out_dir, "cv_prediction_draws.rds"))

# Performance summaries --------------------------------------------------------

# Summarise posterior metric draws
posterior_metric_summary <- function(x) {

  x <- x[is.finite(x)]
  if (length(x) == 0) {

    return(tibble(
        mean = NA_real_,
        median = NA_real_,
        sd = NA_real_,
        q025 = NA_real_,
        q975 = NA_real_
    ))

  }
  
  tibble(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    q975 = quantile(x, 0.975, na.rm = TRUE, names = FALSE)
  )

}

# Label summary columns
rename_summary_cols <- function(dat, prediction = "typical") {

  dat %>%
    rename_with(
      .fn = ~ paste0(prediction, "_", .x),
      .cols = c(mean, median, sd, q025, q975)
    )

}

# Calculate draw-level performance for attrition counts and risk
# Set a floor in case any are extremes (max(r_obs, 1); max(risk_obs, 1 / N))
make_performance_draws <- function(draws, by_fold = FALSE) {

  group_vars <- if (by_fold) c("fold", "draw") else "draw"

  draws %>%
    mutate(
      risk_obs = r_obs / N,
      risk_pred = count_pred / N,
      count_error = count_pred - r_obs,
      abs_count_error = abs(count_error),
      sq_count_error = count_error^2,
      rel_count_error = count_error / pmax(r_obs, 1),
      abs_rel_count_error = abs(rel_count_error),
      risk_error = risk_pred - risk_obs,
      abs_risk_error = abs(risk_error),
      sq_risk_error = risk_error^2,
      rel_risk_error = risk_error / pmax(risk_obs, 1 / N),
      abs_rel_risk_error = abs(rel_risk_error)
    ) %>%
    group_by(across(all_of(group_vars))) %>%
    reframe(
      n = n(),

      mean_obs_count = mean(r_obs),
      mean_pred_count = mean(count_pred),
      total_obs_count = sum(r_obs),
      total_pred_count = sum(count_pred),
      total_count_error = sum(count_error),
      mpe_count = mean(count_error),
      mae_count = mean(abs_count_error),
      rmse_count = sqrt(mean(sq_count_error)),
      mean_rel_error_count = mean(rel_count_error),
      mean_abs_rel_error_count = mean(abs_rel_count_error),
      median_abs_rel_error_count = median(abs_rel_count_error),
      cor_count = safe_cor(r_obs, count_pred),

      # Equal weight per control arm; MPE > 0 means overprediction
      mean_obs_risk = mean(risk_obs),
      mean_pred_risk = mean(risk_pred),
      mpe_risk = mean(risk_error),
      mae_risk = mean(abs_risk_error),
      rmse_risk = sqrt(mean(sq_risk_error)),
      mean_rel_error_risk = mean(rel_risk_error),
      mean_abs_rel_error_risk = mean(abs_rel_risk_error),
      median_abs_rel_error_risk = median(abs_rel_risk_error),
      cor_risk = safe_cor(risk_obs, risk_pred),

      # Sample-size-weighted aggregate risks, distinct from mean arm risks
      aggregate_obs_risk = sum(r_obs) / sum(N),
      aggregate_pred_risk = sum(count_pred) / sum(N),
      aggregate_risk_error = sum(count_error) / sum(N)
    )

}

# Summarise performance across posterior draws among held-out arms
# Intervals describe posterior uncertainty conditional on the observed data
# Performance for a typical held-out trial, i.e. average random effects
make_draw_summary <- function(draws, prediction = "typical") {

  make_performance_draws(draws) %>%
    select(-draw) %>%
    pivot_longer(everything(), names_to = "sum", values_to = "value") %>%
    group_by(sum) %>%
    reframe(posterior_metric_summary(value)) %>%
    rename_summary_cols(prediction)

}

# Generate overall performance summary
cv_sum <- make_draw_summary(
  filter(cv_draws, prediction == "typical"),
  prediction = "typical"
)

write_csv(cv_sum, paste0(z_out_dir, "cv_sum.csv"))

# Fold-specific summaries ----------------------------------------------------

# Summarise draw-level performance by fold
make_fold_draw_summary <- function(draws, prediction = "typical") {

  make_performance_draws(draws, by_fold = TRUE) %>%
    select(-draw) %>%
    pivot_longer(
      cols = -fold,
      names_to = "sum",
      values_to = "value"
    ) %>%
    group_by(fold, sum) %>%
    reframe(posterior_metric_summary(value)) %>%
    rename_summary_cols(prediction)

}

# Generate fold-specific performance summary for typical predictions
cv_sum_fold <- make_fold_draw_summary(
  filter(cv_draws, prediction == "typical"),
  prediction = "typical"
)

# Save fold-specific performance summary
write_csv(cv_sum_fold, paste0(z_out_dir, "cv_sum_folds.csv"))

# Get mean comparison for each fold
cv_avg_folds <- cv_df %>%
  group_by(fold) %>%
  reframe(
    n_trials = n(),
    mean_pred_count = mean(count_typical_mean),
    mean_obs_count = mean(r_obs),
    mean_pred_risk = mean(risk_typical_mean),
    mean_obs_risk = mean(risk_obs)
  )

write_csv(cv_avg_folds, paste0(z_out_dir, "cv_sum_average_folds.csv"))

# Generate count concordance and calibration summaries -------------------------

# Reshape typical trial draws
prediction_ref <- reshape_prediction_draws(filter(
    cv_draws,
    prediction == "typical"
))

# Generate overall concordance summary
count_c_summary <- make_c_summary(
  prediction = "typical",
  pred_obj = prediction_ref
)

# Generate fold-specific concordance summary
count_c_fold_summary <- make_c_fold_summary(
  prediction = "typical",
  pred_obj = prediction_ref
)

# Generate overall calibration summary
count_calibration_summary <- make_calibration_summary(
  prediction = "typical",
  pred_obj = prediction_ref
)

# Generate fold-specific calibration summary
count_calibration_fold_summary <- make_calibration_fold_summary(
  prediction = "typical",
  pred_obj = prediction_ref
)

write_csv(count_c_summary, paste0(z_out_dir, "count_concordance_sum.csv"))
write_csv(count_c_fold_summary, paste0(z_out_dir, "count_concordance_folds.csv"))
write_csv(count_calibration_summary, paste0(z_out_dir, "count_calibration_sum.csv"))
write_csv(count_calibration_fold_summary, paste0(z_out_dir, "count_calibration_folds.csv"))

# Generate count concordance and calibration summaries -------------------------

# Equal-weighting risk concordance
risk_c_summary <- make_c_summary(
  "typical",
  prediction_ref,
  scale = "risk"
)

risk_c_fold_summary <- make_c_fold_summary(
  "typical",
  prediction_ref,
  scale = "risk"
)

# Sample-size-weighted risk concordance
risk_c_size_weighted <- bind_rows(
  make_c_summary("typical", prediction_ref, "risk", size_weighted = TRUE)
)

risk_c_size_weighted_fold <- bind_rows(
  make_c_fold_summary(
    "typical",
    prediction_ref,
    "risk",
    size_weighted = TRUE
  )
)

# Calibration
risk_calibration_summary <- make_calibration_summary(
  "typical",
  prediction_ref,
  scale = "risk"
)
risk_calibration_fold_summary <- make_calibration_fold_summary(
  "typical",
  prediction_ref,
  scale = "risk"
)

write_csv(risk_c_summary, paste0(z_out_dir, "risk_concordance_sum.csv"))
write_csv(risk_c_fold_summary, paste0(z_out_dir, "risk_concordance_folds.csv"))
write_csv(risk_c_size_weighted, paste0(z_out_dir, "risk_concordance_sample_weighted_sum.csv"))
write_csv(risk_c_size_weighted_fold, paste0(z_out_dir, "risk_concordance_sample_weighted_folds.csv"))
write_csv(risk_calibration_summary, paste0(z_out_dir, "risk_calibration_sum.csv"))
write_csv(risk_calibration_fold_summary, paste0(z_out_dir, "risk_calibration_folds.csv"))
