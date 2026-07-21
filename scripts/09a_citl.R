
# Perform calibration-in-the-large for conduct-stage model
# Align participant-marginalised conduct-stage risk to design-stage marginal risk

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Import Sobol outputs
imp_sobol <- readRDS("processed_data/sobol_lst.rds")
imp_sobol_profiles <- readRDS("processed_data/sobol_unnest.rds")

# Import design-stage model and Stan data
imp_design <- readRDS("outputs/design/design_model.rds")
imp_stan <- readRDS("processed_data/stan_data.rds")

# Import AgD
imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Trial quantities
I <- imp_stan$I
r <- imp_stan$r
N <- imp_stan$N
dt <- imp_stan$dt

# Number of points used to invert participant-averaged survival
n_citl_grid <- 4096L

# Check inputs -----------------------------------------------------------------

stopifnot(
  I == nrow(imp_agd),
  I == imp_sobol$I,
  I == nrow(imp_sobol_profiles),
  identical(imp_agd$trial_id, imp_sobol_profiles$trial_id),
  length(r) == I,
  length(N) == I,
  length(dt) == I,
  all(is.finite(r)),
  all(is.finite(N)),
  all(is.finite(dt))
)

# Prepare design-stage posterior draws ----------------------------------------

# Extract Gompertz parameters
design_draws <- as_draws(imp_design)
design_logb <- extract_variable_array(design_draws, variable = "log_b")
design_a <- extract_variable_array(design_draws, variable = "a")

# Convert posterior arrays to matrices
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

logb_mat <- arr_to_mat(design_logb)
a_mat <- arr_to_mat(design_a)

# Number of posterior draws
S <- nrow(logb_mat)

stopifnot(
  ncol(logb_mat) == I,
  all(dim(a_mat) == c(S, I)),
  all(is.finite(logb_mat)),
  all(is.finite(a_mat))
)

# Prepare Sobol participant linear predictors ---------------------------------

J <- imp_sobol$J
Z <- imp_sobol$Z
loghr <- imp_sobol$loghr

# Check predictor names and dimensions
stopifnot(
  setequal(names(loghr), colnames(Z)),
  nrow(Z) == I * J,
  ncol(Z) == length(loghr),
  all(is.finite(Z)),
  all(is.finite(loghr))
)

# Reorder predictors
Z <- Z[, names(loghr), drop = FALSE]

# Get trial-specific participant linear predictors
sobol_lp <- as.vector(Z %*% loghr)
lp_mat <- matrix(sobol_lp, nrow = J, ncol = I)
hr_mat <- exp(lp_mat)

stopifnot(
  all(is.finite(lp_mat)),
  all(is.finite(hr_mat))
)

# Helper functions -------------------------------------------------------------

# Gompertz cumulative hazard
gompertz_cumhaz <- function(b, a, dt) {
  
  out <- numeric(length(b))
  is_exp <- abs(a) < 1e-4
  
  out[is_exp] <- b[is_exp] * dt
  out[!is_exp] <- (
    b[!is_exp] / a[!is_exp]
  ) * expm1(a[!is_exp] * dt)
  
  out
}

# Participant-averaged survival
laplace_survival <- function(h, hr) {
  rowMeans(exp(-outer(h, hr)))
}

# Create interpolation function for CITL shift
make_delta_function <- function(H, hr, n_grid = n_citl_grid, buffer = 6) {
  
  H <- H[is.finite(H) & H > 0]
  
  log_h_grid <- seq(
    log(min(H)) - buffer,
    log(max(H)) + buffer,
    length.out = n_grid
  )
  
  survival_grid <- laplace_survival(exp(log_h_grid), hr)
  marginal_h_grid <- -log(pmax(survival_grid, .Machine$double.xmin))
  
  # Remove duplicated values following numerical underflow
  keep <- is.finite(marginal_h_grid) & !duplicated(marginal_h_grid)
  
  delta_interp <- approxfun(
    marginal_h_grid[keep],
    log_h_grid[keep],
    rule = 2
  )
  
  function(H) {
    delta_interp(H) - log(H)
  }
}

# Backsolve CITL shift ---------------------------------------------------------

target_risk <- matrix(NA_real_, nrow = S, ncol = I)
delta_mat <- matrix(NA_real_, nrow = S, ncol = I)

avg_risk_pre <- matrix(NA_real_, nrow = S, ncol = I)
avg_risk_citl <- matrix(NA_real_, nrow = S, ncol = I)

for (i in seq_len(I)) {
  
  if (i %% 10 == 1 || i == I) {
    message(sprintf("Trial %d/%d", i, I))
  }
  
  # Design-stage cumulative hazard and marginal risk
  H_i <- gompertz_cumhaz(
    b = exp(logb_mat[, i]),
    a = a_mat[, i],
    dt = dt[i]
  )
  
  target_risk[, i] <- -expm1(-H_i)
  
  # Participant-marginalised risk before CITL
  avg_risk_pre[, i] <- 1 - laplace_survival(
    H_i,
    hr_mat[, i]
  )
  
  # Backsolve CITL shift for trial-specific participant population
  delta_fun <- make_delta_function(
    H = H_i,
    hr = hr_mat[, i]
  )
  
  delta_mat[, i] <- delta_fun(H_i)
  
  # Participant-marginalised risk after CITL
  avg_risk_citl[, i] <- 1 - laplace_survival(
    H_i * exp(delta_mat[, i]),
    hr_mat[, i]
  )
}

# Predicted counts
pred_count_design <- sweep(target_risk, 2, N, "*")
pred_count_pre <- sweep(avg_risk_pre, 2, N, "*")
pred_count_citl <- sweep(avg_risk_citl, 2, N, "*")

# Check outputs ----------------------------------------------------------------

citl_resid <- avg_risk_citl - target_risk

stopifnot(
  !anyNA(target_risk),
  !anyNA(delta_mat),
  !anyNA(avg_risk_pre),
  !anyNA(avg_risk_citl),
  max(abs(citl_resid)) < 1e-4
)

summary(as.numeric(citl_resid))

# Save full outputs ------------------------------------------------------------

dir.create("outputs/conduct/citl", recursive = TRUE, showWarnings = FALSE)

citl_out <- list(
  I = I,
  J = J,
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

saveRDS(delta_mat, "outputs/conduct/citl/delta_draws.rds")
saveRDS(citl_out, "outputs/conduct/citl/citl_draws.rds")

# Summarise CITL outputs -------------------------------------------------------

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

# Summarise draw-level metric
summ_draw_metric <- function(x, metric) {
  summ_vec(rowMeans(x, na.rm = TRUE), metric)
}

# Summarise all draws and trials
summ_all_metric <- function(x, metric) {
  summ_vec(as.numeric(x), metric)
}

# Draw-level summaries
citl_draw_summary <- bind_rows(
  summ_draw_metric(target_risk, "design_stage_risk"),
  summ_draw_metric(avg_risk_pre, "pre_citl_risk"),
  summ_draw_metric(avg_risk_citl, "post_citl_risk"),
  summ_draw_metric(avg_risk_pre - target_risk, "pre_minus_design_risk"),
  summ_draw_metric(avg_risk_citl - target_risk, "post_minus_design_risk"),
  summ_draw_metric(delta_mat, "delta_log_rate"),
  summ_draw_metric(pred_count_design, "design_stage_count"),
  summ_draw_metric(pred_count_pre, "pre_citl_count"),
  summ_draw_metric(pred_count_citl, "post_citl_count"),
  summ_draw_metric(pred_count_pre - pred_count_design, "pre_minus_design_count"),
  summ_draw_metric(pred_count_citl - pred_count_design, "post_minus_design_count")
)

# Global summaries
citl_global_summary <- bind_rows(
  summ_all_metric(target_risk, "design_stage_risk"),
  summ_all_metric(avg_risk_pre, "pre_citl_risk"),
  summ_all_metric(avg_risk_citl, "post_citl_risk"),
  summ_all_metric(avg_risk_pre - target_risk, "pre_minus_design_risk"),
  summ_all_metric(avg_risk_citl - target_risk, "post_minus_design_risk"),
  summ_all_metric(delta_mat, "delta_log_rate"),
  summ_all_metric(pred_count_design, "design_stage_count"),
  summ_all_metric(pred_count_pre, "pre_citl_count"),
  summ_all_metric(pred_count_citl, "post_citl_count"),
  summ_all_metric(pred_count_pre - pred_count_design, "pre_minus_design_count"),
  summ_all_metric(pred_count_citl - pred_count_design, "post_minus_design_count")
)

# Save overall summaries
write_csv(citl_draw_summary, "outputs/conduct/citl/citl_draw_summary.csv")
write_csv(citl_global_summary, "outputs/conduct/citl/citl_global_summary.csv")

# Summarise per trial ----------------------------------------------------------

# Summarise posterior draws by trial
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

# Combine trial summaries
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
    arm_summ(avg_risk_pre - target_risk, "risk_pre_minus_design"),
    arm_summ(avg_risk_citl - target_risk, "risk_citl_minus_design"),
    arm_summ(delta_mat, "delta"),
    arm_summ(pred_count_design, "count_design"),
    arm_summ(pred_count_pre, "count_pre"),
    arm_summ(pred_count_citl, "count_citl"),
    arm_summ(pred_count_pre - pred_count_design, "count_pre_minus_design"),
    arm_summ(pred_count_citl - pred_count_design, "count_citl_minus_design")
  )

write_csv(citl_sum, "outputs/conduct/citl/citl_risk_sum.csv")
write_csv(citl_sum, "outputs/conduct/citl/citl_arm_summary.csv")

# Summarise across trial-level posterior medians
summ_arm_medians <- function(data, column, metric) {
  summ_vec(data[[column]], metric)
}

citl_arm_median_summary <- bind_rows(
  summ_arm_medians(citl_sum, "risk_design_median", "arm_median_design_stage_risk"),
  summ_arm_medians(citl_sum, "risk_pre_median", "arm_median_pre_citl_risk"),
  summ_arm_medians(citl_sum, "risk_citl_median", "arm_median_post_citl_risk"),
  summ_arm_medians(citl_sum, "risk_pre_minus_design_median", "arm_median_pre_minus_design_risk"),
  summ_arm_medians(citl_sum, "risk_citl_minus_design_median", "arm_median_post_minus_design_risk"),
  summ_arm_medians(citl_sum, "delta_median", "arm_median_delta_log_rate"),
  summ_arm_medians(citl_sum, "count_design_median", "arm_median_design_stage_count"),
  summ_arm_medians(citl_sum, "count_pre_median", "arm_median_pre_citl_count"),
  summ_arm_medians(citl_sum, "count_citl_median", "arm_median_post_citl_count"),
  summ_arm_medians(citl_sum, "count_pre_minus_design_median", "arm_median_pre_minus_design_count"),
  summ_arm_medians(citl_sum, "count_citl_minus_design_median", "arm_median_post_minus_design_count")
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
  "outputs/conduct/citl/citl_arm_median_summary.csv"
)

write_csv(
  citl_delta_direction,
  "outputs/conduct/citl/citl_delta_direction_summary.csv"
)

# Calibration checks -----------------------------------------------------------

# Summarise residuals
citl_resid_sum <- tibble(
  min = min(citl_resid, na.rm = TRUE),
  q25 = quantile(citl_resid, 0.25, na.rm = TRUE, names = FALSE),
  median = median(citl_resid, na.rm = TRUE),
  mean = mean(citl_resid, na.rm = TRUE),
  q75 = quantile(citl_resid, 0.75, na.rm = TRUE, names = FALSE),
  max = max(citl_resid, na.rm = TRUE),
  n_na = sum(is.na(citl_resid))
)

write_csv(citl_resid_sum, "outputs/conduct/citl/citl_residual_summary.csv")

# Observed-risk diagnostics ----------------------------------------------------

# Compare observed and predicted trial risks
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

# Summarise observed-risk diagnostics
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
  "outputs/conduct/citl/citl_observed_diagnostic_summary.csv"
)

# Plot conduct-stage risk against design-stage target -------------------------

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
      "pre-CITL" = "Pre-CITL",
      "post-CITL" = "Post-CITL"
    ),
    stage = factor(
      stage,
      levels = c("Pre-CITL", "Post-CITL")
    )
  ) %>% 
  ggplot(aes(x = design, y = fit, colour = stage)) +
  geom_abline(linetype = 2, linewidth = 0.5) +
  geom_point(size = 0.5) +
  labs(
    x = "Design-stage marginal attrition risk",
    y = "Participant-marginalised conduct-stage attrition risk",
    colour = "Stage"
  ) +
  theme_classic()

plot_citl

ggsave(
  "outputs/conduct/citl/citl_design_alignment.png",
  plot_citl,
  width = 4,
  height = 4,
  units = "in"
)

# Plot predicted against observed risk ----------------------------------------

plot_obs <- citl_sum %>% 
  select(
    arm,
    trial_id,
    risk_obs,
    risk_design_median,
    risk_pre_median,
    risk_citl_median
  ) %>% 
  pivot_longer(
    cols = c(risk_design_median, risk_pre_median, risk_citl_median),
    names_to = "stage",
    values_to = "risk"
  ) %>% 
  ggplot(aes(x = risk_obs, y = risk, colour = stage)) +
  geom_abline(linetype = 2) +
  geom_point() +
  labs(
    x = "Observed attrition risk",
    y = "Predicted attrition risk"
  ) +
  theme_classic()

plot_obs

ggsave(
  "outputs/conduct/citl/citl_observed_diagnostic.png",
  plot_obs,
  width = 6,
  height = 6,
  units = "in"
)
