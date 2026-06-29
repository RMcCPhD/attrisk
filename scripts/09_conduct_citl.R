
# Perform calibration-in-the-large for the conduct-stage model
# Aligns participant-marginalised conduct-stage risk to the design-stage
# marginal attrition risk after applying pooled participant-level effects
# Finds an intercept shift that aligns observed risk with average participant
# level attrition risk

source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Import Sobol outputs
imp_sobol <- readRDS("processed_data/sobol_lst.rds")

# Import design-stage model and Stan data
imp_design <- readRDS("outputs/design/design_model.rds")
imp_stan <- readRDS("processed_data/stan_data.rds")

# Import AgD data
imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Trial arm quantities
I <- imp_stan$I
r <- imp_stan$r
N <- imp_stan$N
dt <- imp_stan$dt

# Prepare design-stage posterior draws -----------------------------------------

# Extract posterior draws of Gompertz parameters
design_draws <- as_draws(imp_design)
design_logb <- extract_variable_array(design_draws, variable = "log_b")
design_a <- extract_variable_array(design_draws, variable = "a")

# Convert posterior draw arrays to matrices (draws x trial arms)
arr_to_mat <- function(arr) {
  d <- dim(arr)
  if (length(d) == 3) {
    matrix(arr, nrow = d[1] * d[2], ncol = d[3])
  } else if (length(d) == 2) {
    as.matrix(arr)
  } else {
    stop("Unexpected posterior draw dimension.")
  }
}

logb_mat <- arr_to_mat(design_logb)
a_mat <- arr_to_mat(design_a)

# Prepare Sobol participant linear predictors ----------------------------------

J <- imp_sobol$J
Z <- imp_sobol$Z
loghr <- imp_sobol$loghr

# Checks names match
setdiff(names(loghr), colnames(Z))
setdiff(colnames(Z), names(loghr))

# Reorder Z block explicitly to names
Z <- Z[, names(loghr), drop = FALSE]

# One participant-level linear predictor per Sobol profile
sobol_lp <- as.vector(Z %*% loghr)

# Split linear predictors into a vector of J participants per trial arm
lp_lst <- split(sobol_lp, rep(seq_len(I), each = J))

# Gompertz CDF -----------------------------------------------------------------

# If shape is close to zero, decays to exponential
gompertz_prob <- function(b, a, dt, lp = 0) {
  if (abs(a) < 1e-4) {
    p <- 1 - exp(-b * exp(lp) * dt)
  } else {
    p <- 1 - exp(-(b / a) * (exp(a * dt) - 1) * exp(lp))
  }
  pmin(pmax(p, 0), 1)
}

# Root-finding helpers ---------------------------------------------------------

uniroot_bracket <- function(f, lo = -30, hi = 30, step = 5, max_iter = 10) {
  
  flo <- f(lo)
  fhi <- f(hi)
  k <- 0L
  
  if (!is.finite(flo) || !is.finite(fhi)) return(NULL)
  if (abs(flo) < 1e-12) return(lo)
  if (abs(fhi) < 1e-12) return(hi)
  
  while (sign(flo) == sign(fhi) && k < max_iter) {
    lo <- lo - step
    hi <- hi + step
    flo <- f(lo)
    fhi <- f(hi)
    k <- k + 1L
    
    if (!is.finite(flo) || !is.finite(fhi)) return(NULL)
    if (abs(flo) < 1e-12) return(lo)
    if (abs(fhi) < 1e-12) return(hi)
  }
  
  if (sign(flo) == sign(fhi)) return(NULL)
  
  uniroot(f, c(lo, hi))$root
  
}

# Find the intercept shift that aligns the participant-marginalised risk with
# the design-stage marginal risk.
find_delta <- function(log_b, a, dt, lp, target_risk) {
  
  target_risk <- pmin(pmax(target_risk, 1e-10), 1 - 1e-10)
  
  f <- function(delta) {
    p <- gompertz_prob(
      b = exp(log_b + delta),
      a = a,
      dt = dt,
      lp = lp
    )
    mean(p) - target_risk
  }
  
  out <- uniroot_bracket(f, lo = -30, hi = 30, step = 5)
  
  if (is.null(out)) {
    obj <- function(delta) (f(delta))^2
    optimise(obj, interval = c(-40, 40))$minimum
  } else {
    out
  }
  
}

# Backsolve to find delta (CITL intercept shift) -------------------------------

S <- nrow(logb_mat)

target_risk <- matrix(NA_real_, nrow = S, ncol = I)
delta_mat <- matrix(NA_real_, nrow = S, ncol = I)

for (i in seq_len(I)) {
  
  lp_i <- lp_lst[[i]]
  dt_i <- dt[i]
  
  if (i %% 10 == 1 || i == I) message(sprintf("Arm %d/%d", i, I))
  
  for (s in seq_len(S)) {
    
    # Design-stage marginal risk before applying participant-level effects
    target_risk[s, i] <- gompertz_prob(
      b = exp(logb_mat[s, i]),
      a = a_mat[s, i],
      dt = dt_i,
      lp = 0
    )
    
    # CITL shift recentres the participant-marginalised risk to this target
    delta_mat[s, i] <- find_delta(
      log_b = logb_mat[s, i],
      a = a_mat[s, i],
      dt = dt_i,
      lp = lp_i,
      target_risk = target_risk[s, i]
    )
    
  }
  
}

mean(is.na(delta_mat))

# Compute pre- and post-CITL risks --------------------------------------------

avg_risk_pre <- matrix(NA_real_, nrow = S, ncol = I)
pred_count_pre <- matrix(NA_real_, nrow = S, ncol = I)

avg_risk_citl <- matrix(NA_real_, nrow = S, ncol = I)
pred_count_citl <- matrix(NA_real_, nrow = S, ncol = I)

for (i in seq_len(I)) {
  
  lp_i <- lp_lst[[i]]
  dt_i <- dt[i]
  
  if (i %% 10 == 1 || i == I) message(sprintf("Arm %d/%d", i, I))
  
  for (s in seq_len(S)) {
    
    p_pre <- gompertz_prob(
      b = exp(logb_mat[s, i]),
      a = a_mat[s, i],
      dt = dt_i,
      lp = lp_i
    )
    
    p_citl <- gompertz_prob(
      b = exp(logb_mat[s, i] + delta_mat[s, i]),
      a = a_mat[s, i],
      dt = dt_i,
      lp = lp_i
    )
    
    avg_risk_pre[s, i] <- mean(p_pre)
    avg_risk_citl[s, i] <- mean(p_citl)
    
    pred_count_pre[s, i] <- N[i] * avg_risk_pre[s, i]
    pred_count_citl[s, i] <- N[i] * avg_risk_citl[s, i]
    
  }
  
}

# Save full outputs ------------------------------------------------------------

dir.create("outputs/conduct/citl", recursive = TRUE, showWarnings = FALSE)

citl_out <- list(
  I = I,
  J = J,
  S = S,
  target_risk = target_risk,
  delta = delta_mat,
  avg_risk_pre = avg_risk_pre,
  avg_risk_citl = avg_risk_citl,
  pred_count_pre = pred_count_pre,
  pred_count_citl = pred_count_citl
)

saveRDS(delta_mat, "outputs/conduct/citl/delta_draws.rds")
saveRDS(citl_out, "outputs/conduct/citl/citl_draws.rds")

# Summarise CITL outputs -------------------------------------------------------

summ_vec <- function(x, metric) {
  tibble(
    metric = metric,
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE),
    q975 = quantile(x, 0.975, na.rm = TRUE)
  )
}

summ_draw_metric <- function(x, metric) {
  summ_vec(rowMeans(x, na.rm = TRUE), metric)
}

summ_all_metric <- function(x, metric) {
  summ_vec(as.numeric(x), metric)
}

count_design <- sweep(target_risk, 2, N, "*")
count_pre <- pred_count_pre
count_citl <- pred_count_citl

citl_out$pred_count_design <- count_design
saveRDS(citl_out, "outputs/conduct/citl/citl_draws.rds")

# Draw-level summaries: each posterior draw is first averaged across AgD arms.
citl_draw_summary <- bind_rows(
  summ_draw_metric(target_risk, "design_stage_risk"),
  summ_draw_metric(avg_risk_pre, "pre_citl_risk"),
  summ_draw_metric(avg_risk_citl, "post_citl_risk"),
  summ_draw_metric(avg_risk_pre - target_risk, "pre_minus_design_risk"),
  summ_draw_metric(avg_risk_citl - target_risk, "post_minus_design_risk"),
  summ_draw_metric(delta_mat, "delta_log_rate"),
  summ_draw_metric(count_design, "design_stage_count"),
  summ_draw_metric(count_pre, "pre_citl_count"),
  summ_draw_metric(count_citl, "post_citl_count"),
  summ_draw_metric(count_pre - count_design, "pre_minus_design_count"),
  summ_draw_metric(count_citl - count_design, "post_minus_design_count")
)

# Global summaries: all posterior draws and arms are pooled before summarising.
citl_global_summary <- bind_rows(
  summ_all_metric(target_risk, "design_stage_risk"),
  summ_all_metric(avg_risk_pre, "pre_citl_risk"),
  summ_all_metric(avg_risk_citl, "post_citl_risk"),
  summ_all_metric(avg_risk_pre - target_risk, "pre_minus_design_risk"),
  summ_all_metric(avg_risk_citl - target_risk, "post_minus_design_risk"),
  summ_all_metric(delta_mat, "delta_log_rate"),
  summ_all_metric(count_design, "design_stage_count"),
  summ_all_metric(count_pre, "pre_citl_count"),
  summ_all_metric(count_citl, "post_citl_count"),
  summ_all_metric(count_pre - count_design, "pre_minus_design_count"),
  summ_all_metric(count_citl - count_design, "post_minus_design_count")
)

write_csv(citl_draw_summary, "outputs/conduct/citl/citl_draw_summary.csv")
write_csv(citl_global_summary, "outputs/conduct/citl/citl_global_summary.csv")

# Summarise per trial arm ------------------------------------------------------

arm_summ <- function(x, prefix) {
  tibble(
    "{prefix}_mean" := apply(x, 2, mean, na.rm = TRUE),
    "{prefix}_median" := apply(x, 2, median, na.rm = TRUE),
    "{prefix}_sd" := apply(x, 2, sd, na.rm = TRUE),
    "{prefix}_q025" := apply(x, 2, \(z) quantile(z, 0.025, na.rm = TRUE)),
    "{prefix}_q975" := apply(x, 2, \(z) quantile(z, 0.975, na.rm = TRUE))
  )
}

citl_sum <- tibble(
  arm = seq_len(I),
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
    arm_summ(count_design, "count_design"),
    arm_summ(count_pre, "count_pre"),
    arm_summ(count_citl, "count_citl"),
    arm_summ(count_pre - count_design, "count_pre_minus_design"),
    arm_summ(count_citl - count_design, "count_citl_minus_design")
  )

write_csv(citl_sum, "outputs/conduct/citl/citl_risk_sum.csv")
write_csv(citl_sum, "outputs/conduct/citl/citl_arm_summary.csv")

# Empirical summaries across arm-level posterior medians, useful for text.
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

# Should be very close to zero
citl_resid <- avg_risk_citl - target_risk

summary(as.numeric(citl_resid))

citl_resid_sum <- tibble(
  min = min(citl_resid, na.rm = TRUE),
  q25 = quantile(citl_resid, 0.25, na.rm = TRUE),
  median = median(citl_resid, na.rm = TRUE),
  mean = mean(citl_resid, na.rm = TRUE),
  q75 = quantile(citl_resid, 0.75, na.rm = TRUE),
  max = max(citl_resid, na.rm = TRUE),
  n_na = sum(is.na(citl_resid))
)

write_csv(citl_resid_sum, "outputs/conduct/citl/citl_residual_summary.csv")

# Observed-risk diagnostics ----------------------------------------------------

observed_diag_dat <- citl_sum %>%
  transmute(
    arm,
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
  "outputs/conduct/citl/citl_observed_diagnostic_summary.csv"
)

# Plot: conduct-stage risk against design-stage target -------------------------

plot_dat <- tibble(
  arm = rep(seq_len(I), 2),
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

# Diagnostic plot against observed risk ----------------------------------------

plot_obs <- citl_sum %>% 
  select(arm, risk_obs, risk_design_median, risk_pre_median, risk_citl_median) %>% 
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
