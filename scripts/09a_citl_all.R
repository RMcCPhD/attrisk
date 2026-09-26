
# Perform CITL across all completed AgD datasets
# Use the pooled participant-profile distribution as the primary input

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Import Sobol outputs for all AgD imputations
imp_sobol_mi <- readRDS("outputs/conduct_stage/profiles_imputation.rds")

# Import design-stage model and Stan data
imp_design <- readRDS("outputs/design_stage/model/design_model.rds")
imp_stan <- readRDS("processed_data/stan_data.rds")

# Import AgD
imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Output directory
z_out <- "outputs/conduct_stage/citl"
dir.create(z_out, recursive = TRUE, showWarnings = FALSE)

# Trial quantities
I <- imp_stan$I
r <- imp_stan$r
N <- imp_stan$N
dt <- imp_stan$dt

# Number of points used to invert participant-averaged survival
n_citl_grid <- 4096L

# Number of Sobol profiles per arm in each completed dataset
n_sobol_per_imputation <- 64L

# Helper functions -------------------------------------------------------------

# Convert posterior array to matrix
arr_to_mat <- function(arr) {
  d <- dim(arr)
  if (length(d) == 3) {
    matrix(arr, nrow = d[1] * d[2], ncol = d[3])
  } else if (length(d) == 2) {
    as.matrix(arr)
  } else {
    stop("Unexpected posterior draw dimension")
  }
}

# Gompertz cumulative hazard
gompertz_cumhaz <- function(b, a, dt) {
  out <- numeric(length(b))
  is_exp <- abs(a) < 1e-4
  out[is_exp] <- b[is_exp] * dt
  out[!is_exp] <- (b[!is_exp] / a[!is_exp]) * expm1(a[!is_exp] * dt)
  out
}

# Participant-averaged survival
mean_survival <- function(h, hr) {rowMeans(exp(-outer(h, hr)))}

# Create interpolation function for CITL shift
make_delta_function <- function(H, hr, n_grid = n_citl_grid, buffer = 6) {
  H <- H[is.finite(H) & H > 0]
  log_h_grid <- seq(
    log(min(H)) - buffer,
    log(max(H)) + buffer,
    length.out = n_grid
  )
  survival_grid <- mean_survival(exp(log_h_grid), hr)
  marginal_h_grid <- -log(pmax(survival_grid, .Machine$double.xmin))
  keep <- is.finite(marginal_h_grid) & !duplicated(marginal_h_grid)
  delta_interp <- approxfun(
    marginal_h_grid[keep],
    log_h_grid[keep],
    rule = 2
  )
  function(H) {delta_interp(H) - log(H)}
}

# Summarise vector
summ_vec <- function(x, metric) {
  tibble(
    metric = metric,
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    q975 = quantile(x, 0.975, na.rm = TRUE, names = FALSE)
  )
}

# Summarise posterior matrix by arm
arm_summ <- function(x, prefix) {
  tibble(
    "{prefix}_mean" := apply(x, 2, mean, na.rm = TRUE),
    "{prefix}_median" := apply(x, 2, median, na.rm = TRUE),
    "{prefix}_sd" := apply(x, 2, sd, na.rm = TRUE),
    "{prefix}_q025" := apply(
      x,
      2,
      \(z) quantile(z, 0.025, na.rm = TRUE, names = FALSE)
    ),
    "{prefix}_q975" := apply(
      x,
      2,
      \(z) quantile(z, 0.975, na.rm = TRUE, names = FALSE)
    )
  )
}

# Convert Sobol input to profile linear-predictor matrix
get_lp_matrix <- function(sobol) {
  Z <- sobol$Z[, names(sobol$loghr), drop = FALSE]
  lp <- as.vector(Z %*% sobol$loghr)
  out <- matrix(
    lp,
    nrow = sobol$J,
    ncol = sobol$I
  )
  stopifnot(
    ncol(out) == I,
    all(is.finite(out))
  )
  out
}

# Prepare design-stage posterior draws -----------------------------------------

design_draws <- as_draws(imp_design)
design_logb <- extract_variable_array(design_draws, variable = "log_b")
design_a <- extract_variable_array(design_draws, variable = "a")
logb_mat <- arr_to_mat(design_logb)
a_mat <- arr_to_mat(design_a)
S <- nrow(logb_mat)

stopifnot(
  ncol(logb_mat) == I,
  all(dim(a_mat) == c(S, I)),
  all(is.finite(logb_mat)),
  all(is.finite(a_mat))
)

# Calculate design-stage cumulative hazards once
H_mat <- vapply(
  seq_len(I),
  function(i) {
    gompertz_cumhaz(
      b = exp(logb_mat[, i]),
      a = a_mat[, i],
      dt = dt[i]
    )
  },
  numeric(S)
)

target_risk <- -expm1(-H_mat)
pred_count_design <- sweep(target_risk, 2, N, "*")

stopifnot(
  all(dim(H_mat) == c(S, I)),
  all(is.finite(H_mat)),
  all(H_mat > 0),
  all(is.finite(target_risk))
)

# Perform CITL for one profile set ---------------------------------------------

run_citl_profiles <- function(lp_mat, label) {

  J <- nrow(lp_mat)
  hr_mat <- exp(lp_mat)

  stopifnot(
    ncol(lp_mat) == I,
    all(is.finite(lp_mat)),
    all(is.finite(hr_mat))
  )

  delta_mat <- matrix(NA_real_, nrow = S, ncol = I)
  avg_risk_pre <- matrix(NA_real_, nrow = S, ncol = I)
  avg_risk_citl <- matrix(NA_real_, nrow = S, ncol = I)

  for (i in seq_len(I)) {

    if (i %% 25 == 1 || i == I) {

      message(sprintf("%s: arm %d/%d", label, i, I))

    }

    H_i <- H_mat[, i]
    hr_i <- hr_mat[, i]

    avg_risk_pre[, i] <- 1 - mean_survival(H_i, hr_i)

    delta_fun <- make_delta_function(
      H = H_i,
      hr = hr_i
    )

    delta_mat[, i] <- delta_fun(H_i)

    avg_risk_citl[, i] <- 1 - mean_survival(
      H_i * exp(delta_mat[, i]),
      hr_i
    )

  }

  pred_count_pre <- sweep(avg_risk_pre, 2, N, "*")
  pred_count_citl <- sweep(avg_risk_citl, 2, N, "*")
  citl_resid <- avg_risk_citl - target_risk

  stopifnot(
    !anyNA(delta_mat),
    !anyNA(avg_risk_pre),
    !anyNA(avg_risk_citl),
    max(abs(citl_resid)) < 1e-4
  )

  # Profile summaries by arm
  profile_summary <- tibble(
    arm = seq_len(I),
    trial_id = imp_agd$trial_id,
    n_profiles = J,
    lp_mean = apply(lp_mat, 2, mean),
    lp_median = apply(lp_mat, 2, median),
    lp_sd = apply(lp_mat, 2, sd),
    lp_q025 = apply(lp_mat, 2, quantile, probs = 0.025, names = FALSE),
    lp_q975 = apply(lp_mat, 2, quantile, probs = 0.975, names = FALSE),
    hr_multiplier_mean = apply(hr_mat, 2, mean),
    hr_multiplier_median = apply(hr_mat, 2, median),
    hr_multiplier_sd = apply(hr_mat, 2, sd),
    hr_multiplier_q025 = apply(
      hr_mat,
      2,
      quantile,
      probs = 0.025,
      names = FALSE
    ),
    hr_multiplier_q975 = apply(
      hr_mat,
      2,
      quantile,
      probs = 0.975,
      names = FALSE
    )
  )

  # Posterior summaries by arm
  arm_summary <- tibble(
    arm = seq_len(I),
    trial_id = imp_agd$trial_id,
    N = N,
    r_obs = r,
    risk_obs = r / N
  ) %>%
    left_join(profile_summary, by = c("arm", "trial_id")) %>%
    bind_cols(
      arm_summ(target_risk, "risk_design"),
      arm_summ(avg_risk_pre, "risk_pre"),
      arm_summ(avg_risk_citl, "risk_citl"),
      arm_summ(delta_mat, "delta"),
      arm_summ(pred_count_design, "count_design"),
      arm_summ(pred_count_pre, "count_pre"),
      arm_summ(pred_count_citl, "count_citl")
    )

  # Overall summaries
  global_summary <- bind_rows(
    summ_vec(as.numeric(lp_mat), "participant_linear_predictor"),
    summ_vec(as.numeric(hr_mat), "participant_hazard_multiplier"),
    summ_vec(as.numeric(target_risk), "design_stage_risk"),
    summ_vec(as.numeric(avg_risk_pre), "pre_citl_risk"),
    summ_vec(as.numeric(avg_risk_citl), "post_citl_risk"),
    summ_vec(as.numeric(delta_mat), "delta_log_rate"),
    summ_vec(as.numeric(pred_count_design), "design_stage_count"),
    summ_vec(as.numeric(pred_count_pre), "pre_citl_count"),
    summ_vec(as.numeric(pred_count_citl), "post_citl_count")
  )

  out <- list(
    label = label,
    J = J,
    arm_summary = arm_summary,
    global_summary = global_summary,
    residual_summary = tibble(
      min = min(citl_resid),
      median = median(citl_resid),
      mean = mean(citl_resid),
      max = max(citl_resid),
      max_abs = max(abs(citl_resid))
    )
  )

  out$draws <- list(
    lp = lp_mat,
    hr_multiplier = hr_mat,
    target_risk = target_risk,
    delta = delta_mat,
    avg_risk_pre = avg_risk_pre,
    avg_risk_citl = avg_risk_citl,
    pred_count_design = pred_count_design,
    pred_count_pre = pred_count_pre,
    pred_count_citl = pred_count_citl
  )

  out

}

# Prepare profile sets across imputations --------------------------------------

profile_inputs <- imp_sobol_mi %>%
  arrange(.imp) %>%
  transmute(
    .imp,
    lp_mat = map(sobol_lst, get_lp_matrix)
  )

# The large Sobol design matrices are no longer needed after calculating LPs
rm(imp_sobol_mi)

# Number of integration profiles per arm in each completed dataset
J_by_imputation <- map_int(profile_inputs$lp_mat, nrow)
M <- nrow(profile_inputs)

stopifnot(
  length(J_by_imputation) == M,
  all(J_by_imputation == n_sobol_per_imputation)
)

# Equally weighted mixture of profiles from all completed datasets
pooled_lp <- do.call(rbind, profile_inputs$lp_mat)

stopifnot(
  nrow(pooled_lp) == sum(J_by_imputation),
  ncol(pooled_lp) == I,
  all(is.finite(pooled_lp))
)

# Perform CITL using profiles pooled across imputations ------------------------

citl_pooled <- run_citl_profiles(
  lp_mat = pooled_lp,
  label = sprintf(
    "Pooled imputations (%d profiles per imputation)",
    n_sobol_per_imputation
  )
)

write_csv(
  citl_pooled$arm_summary,
  file.path(z_out, "citl_arm_summary_pooled_profiles.csv")
)

write_csv(
  citl_pooled$global_summary,
  file.path(z_out, "citl_global_summary_pooled_profiles.csv")
)

# Promote pooled-imputation CITL to primary pipeline -------------------------

# Extract CITL draws estimated from profiles pooled across imputations
# Use pooled-profile shifts directly
pooled_draws <- citl_pooled$draws

lp_mat <- pooled_draws$lp
hr_mat <- pooled_draws$hr_multiplier
sobol_lp <- as.vector(lp_mat)

target_risk <- pooled_draws$target_risk
delta_mat <- pooled_draws$delta
avg_risk_pre <- pooled_draws$avg_risk_pre
avg_risk_citl <- pooled_draws$avg_risk_citl

pred_count_design <- pooled_draws$pred_count_design
pred_count_pre <- pooled_draws$pred_count_pre
pred_count_citl <- pooled_draws$pred_count_citl

citl_resid <- avg_risk_citl - target_risk

# Number of profiles per arm after pooling across completed datasets
J <- nrow(lp_mat)

stopifnot(
  length(J_by_imputation) == M,
  all(J_by_imputation == n_sobol_per_imputation),
  J == sum(J_by_imputation),
  ncol(lp_mat) == I,
  nrow(delta_mat) == S,
  ncol(delta_mat) == I,
  all(dim(target_risk) == c(S, I)),
  all(dim(avg_risk_pre) == c(S, I)),
  all(dim(avg_risk_citl) == c(S, I)),
  all(dim(pred_count_design) == c(S, I)),
  all(dim(pred_count_pre) == c(S, I)),
  all(dim(pred_count_citl) == c(S, I)),
  all(is.finite(delta_mat)),
  all(is.finite(avg_risk_pre)),
  all(is.finite(avg_risk_citl)),
  max(abs(citl_resid)) < 1e-4
)

summary(as.numeric(citl_resid))

# Save primary pooled-profile draw outputs -------------------------------------

citl_out <- list(
  I = I,
  J = J,
  J_by_imputation = J_by_imputation,
  n_sobol_per_imputation = n_sobol_per_imputation,
  M = M,
  imputation = profile_inputs$.imp,
  S = S,
  trial_id = imp_agd$trial_id,
  sobol_lp = sobol_lp,
  n_citl_grid = n_citl_grid,
  target_risk = target_risk,
  delta = delta_mat,
  avg_risk_pre = avg_risk_pre,
  avg_risk_citl = avg_risk_citl,
  pred_count_design = pred_count_design,
  pred_count_pre = pred_count_pre,
  pred_count_citl = pred_count_citl
)

saveRDS(delta_mat, file.path(z_out, "delta_draws.rds"))
saveRDS(citl_out, file.path(z_out, "citl_draws.rds"))

# Summarise primary pooled-profile CITL outputs --------------------------------

# Summarise an arm-averaged value for every posterior draw
summ_draw_metric <- function(x, metric) {
  summ_vec(rowMeans(x, na.rm = TRUE), metric)
}

# Summarise all posterior draws and control arms together
summ_all_metric <- function(x, metric) {
  summ_vec(as.numeric(x), metric)
}

# Draw-level summaries
citl_draw_summary <- bind_rows(
  summ_draw_metric(target_risk, "design_stage_risk"),
  summ_draw_metric(avg_risk_pre, "pre_citl_risk"),
  summ_draw_metric(avg_risk_citl, "post_citl_risk"),
  summ_draw_metric(
    avg_risk_pre - target_risk,
    "pre_minus_design_risk"
  ),
  summ_draw_metric(
    avg_risk_citl - target_risk,
    "post_minus_design_risk"
  ),
  summ_draw_metric(delta_mat, "delta_log_rate"),
  summ_draw_metric(pred_count_design, "design_stage_count"),
  summ_draw_metric(pred_count_pre, "pre_citl_count"),
  summ_draw_metric(pred_count_citl, "post_citl_count"),
  summ_draw_metric(
    pred_count_pre - pred_count_design,
    "pre_minus_design_count"
  ),
  summ_draw_metric(
    pred_count_citl - pred_count_design,
    "post_minus_design_count"
  )
)

# Global summaries
citl_global_summary <- bind_rows(
  summ_all_metric(lp_mat, "participant_linear_predictor"),
  summ_all_metric(hr_mat, "participant_hazard_multiplier"),
  summ_all_metric(target_risk, "design_stage_risk"),
  summ_all_metric(avg_risk_pre, "pre_citl_risk"),
  summ_all_metric(avg_risk_citl, "post_citl_risk"),
  summ_all_metric(
    avg_risk_pre - target_risk,
    "pre_minus_design_risk"
  ),
  summ_all_metric(
    avg_risk_citl - target_risk,
    "post_minus_design_risk"
  ),
  summ_all_metric(delta_mat, "delta_log_rate"),
  summ_all_metric(pred_count_design, "design_stage_count"),
  summ_all_metric(pred_count_pre, "pre_citl_count"),
  summ_all_metric(pred_count_citl, "post_citl_count"),
  summ_all_metric(
    pred_count_pre - pred_count_design,
    "pre_minus_design_count"
  ),
  summ_all_metric(
    pred_count_citl - pred_count_design,
    "post_minus_design_count"
  )
)

write_csv(
  citl_draw_summary,
  file.path(z_out, "citl_draw_summary.csv")
)

write_csv(
  citl_global_summary,
  file.path(z_out, "citl_global_summary.csv")
)

# Summarise primary pooled-profile outputs by control arm --------------------

citl_sum <- tibble(
  arm = seq_len(I),
  trial_id = imp_agd$trial_id,
  N = N,
  r_obs = r,
  risk_obs = r / N
) %>%
  bind_cols(
    arm_summ(target_risk, "risk_design"),
    arm_summ(avg_risk_pre, "risk_pre"),
    arm_summ(avg_risk_citl, "risk_citl"),
    arm_summ(
      avg_risk_pre - target_risk,
      "risk_pre_minus_design"
    ),
    arm_summ(
      avg_risk_citl - target_risk,
      "risk_citl_minus_design"
    ),
    arm_summ(delta_mat, "delta"),
    arm_summ(pred_count_design, "count_design"),
    arm_summ(pred_count_pre, "count_pre"),
    arm_summ(pred_count_citl, "count_citl"),
    arm_summ(
      pred_count_pre - pred_count_design,
      "count_pre_minus_design"
    ),
    arm_summ(
      pred_count_citl - pred_count_design,
      "count_citl_minus_design"
    )
  )

write_csv(
  citl_sum,
  file.path(z_out, "citl_risk_sum.csv")
)

write_csv(
  citl_sum,
  file.path(z_out, "citl_arm_summary.csv")
)

# Summarise across control-arm posterior medians
summ_arm_medians <- function(data, column, metric) {
  summ_vec(data[[column]], metric)
}

citl_arm_median_summary <- bind_rows(
  summ_arm_medians(
    citl_sum,
    "risk_design_median",
    "arm_median_design_stage_risk"
  ),
  summ_arm_medians(
    citl_sum,
    "risk_pre_median",
    "arm_median_pre_citl_risk"
  ),
  summ_arm_medians(
    citl_sum,
    "risk_citl_median",
    "arm_median_post_citl_risk"
  ),
  summ_arm_medians(
    citl_sum,
    "risk_pre_minus_design_median",
    "arm_median_pre_minus_design_risk"
  ),
  summ_arm_medians(
    citl_sum,
    "risk_citl_minus_design_median",
    "arm_median_post_minus_design_risk"
  ),
  summ_arm_medians(
    citl_sum,
    "delta_median",
    "arm_median_delta_log_rate"
  ),
  summ_arm_medians(
    citl_sum,
    "count_design_median",
    "arm_median_design_stage_count"
  ),
  summ_arm_medians(
    citl_sum,
    "count_pre_median",
    "arm_median_pre_citl_count"
  ),
  summ_arm_medians(
    citl_sum,
    "count_citl_median",
    "arm_median_post_citl_count"
  ),
  summ_arm_medians(
    citl_sum,
    "count_pre_minus_design_median",
    "arm_median_pre_minus_design_count"
  ),
  summ_arm_medians(
    citl_sum,
    "count_citl_minus_design_median",
    "arm_median_post_minus_design_count"
  )
)

# Summarise direction of CITL shifts
citl_delta_direction <- tibble(
  metric = c(
    "n_arms",
    "n_arms_delta_median_below_zero",
    "prop_arms_delta_median_below_zero"
  ),
  value = c(
    I,
    sum(citl_sum$delta_median < 0, na.rm = TRUE),
    mean(citl_sum$delta_median < 0, na.rm = TRUE)
  )
)

write_csv(
  citl_arm_median_summary,
  file.path(z_out, "citl_arm_median_summary.csv")
)

write_csv(
  citl_delta_direction,
  file.path(z_out, "citl_delta_direction_summary.csv")
)

# Calibration checks for primary pooled-profile CITL -------------------------

citl_resid_sum <- tibble(
  min = min(citl_resid, na.rm = TRUE),
  q25 = quantile(citl_resid, 0.25, na.rm = TRUE, names = FALSE),
  median = median(citl_resid, na.rm = TRUE),
  mean = mean(citl_resid, na.rm = TRUE),
  q75 = quantile(citl_resid, 0.75, na.rm = TRUE, names = FALSE),
  max = max(citl_resid, na.rm = TRUE),
  max_abs = max(abs(citl_resid), na.rm = TRUE),
  n_na = sum(is.na(citl_resid))
)

write_csv(
  citl_resid_sum,
  file.path(z_out, "citl_residual_summary.csv")
)

# Observed-risk diagnostics --------------------------------------------------

observed_diag_dat <- citl_sum %>%
  transmute(
    arm,
    trial_id,
    risk_obs,
    `Design stage` = risk_design_median,
    `Pre-CITL` = risk_pre_median,
    `Post-CITL` = risk_citl_median
  ) %>%
  pivot_longer(
    cols = c(`Design stage`, `Pre-CITL`, `Post-CITL`),
    names_to = "stage",
    values_to = "risk_pred"
  )

observed_diag_summary <- observed_diag_dat %>%
  group_by(stage) %>%
  summarise(
    mean_observed_risk = mean(risk_obs, na.rm = TRUE),
    mean_predicted_risk = mean(risk_pred, na.rm = TRUE),
    mean_error = mean(risk_pred - risk_obs, na.rm = TRUE),
    mae = mean(abs(risk_pred - risk_obs), na.rm = TRUE),
    rmse = sqrt(mean((risk_pred - risk_obs)^2, na.rm = TRUE)),
    cor_obs_pred = cor(risk_obs, risk_pred, use = "complete.obs"),
    .groups = "drop"
  )

write_csv(
  observed_diag_summary,
  file.path(z_out, "citl_observed_diagnostic_summary.csv")
)

# Plot conduct-stage risk against design-stage target ------------------------

plot_dat <- tibble(
  arm = rep(seq_len(I), 2),
  trial_id = rep(imp_agd$trial_id, 2),
  design = rep(apply(target_risk, 2, median, na.rm = TRUE), 2),
  fit = c(
    apply(avg_risk_pre, 2, median, na.rm = TRUE),
    apply(avg_risk_citl, 2, median, na.rm = TRUE)
  ),
  stage = rep(c("pre-CITL", "post-CITL"), each = I)
)

plot_citl <- plot_dat %>%
  mutate(
    stage = recode(
      stage,
      "pre-CITL" = "Before",
      "post-CITL" = "After"
    ),
    stage = factor(
      stage,
      levels = c("Before", "After")
    )
  ) %>%
  ggplot(aes(x = design, y = fit, colour = stage)) +
  geom_abline(linetype = 2, linewidth = 0.5) +
  geom_point(size = 0.5) +
  labs(
    x = "Design-stage predicted attrition risk",
    y = "Conduct-stage predicted attrition risk",
    colour = "Alignment stage"
  ) +
  theme_classic()

plot_citl

ggsave(
  file.path(z_out, "citl_calibration_plot.png"),
  plot_citl,
  width = 4,
  height = 4,
  units = "in"
)
