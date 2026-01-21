
# Find delta (i.e. intercept shift) for CITL
# Pull in posterior draws from conduct model, Sobols and pooled HRs
# Backsolve with uniroot to find delta that aligns observed risk in trials
# with average risk among participants

source("scripts/00_config.R")
library(tidyverse)
library(posterior)

# Import pre-CITL model, get draws
a_imp_mdl <- readRDS("mars/conduct_model/outputs/conduct_fit.rds")
a_draws <- as_draws(a_imp_mdl)

# Import pre-CITL dataset
a_imp_standf <- readRDS("mars/design_model/latest_stan_design.rds")

# Import sobol outputs
a_imp_sobol <- readRDS("processed_data/sobol_lst.rds")

# Get posterior draws of log-rate, rate and shape
b_logb <- extract_variable_array(a_draws, variable = "log_b")
b_a <- extract_variable_array(a_draws, variable = "a")

# Reshape arrays into matrices of draws x trials
arr_to_mat <- function(arr) {
  iter <- dim(arr)[1]
  chain <- dim(arr)[2]
  trial <- dim(arr)[3]
  matrix(arr, nrow = iter * chain, ncol = trial)
}

b_logb_mat <- arr_to_mat(b_logb)
b_a_mat <- arr_to_mat(b_a)

# Pull out trials, counts (i.e. target), arm size and time
I <- a_imp_standf$I
r <- a_imp_standf$r
N <- a_imp_standf$N
dt <- a_imp_standf$dt

# Build arm-level participant linear predictors (log-HR), split by arm
b_sobol_lp <- as.vector(a_imp_sobol$Z %*% a_imp_sobol$loghr)
b_lp <- split(b_sobol_lp, rep(seq_len(I), each = 64L))

# Build time list
b_dt_lst <- lapply(seq_len(I), \(i) rep(dt[i], length(b_lp[[i]])))

# Compute probability using Gompertz CDF, limit to exp if shape close to 0
gompertz_prob <- function(b, a, dt, lp) {
  if (abs(a) < 1e-4) {
    1 - exp(-b * exp(lp) * dt)
  } else {
    1 - exp(-(b / a) * (exp(a * dt) - 1) * exp(lp))
  }
}

# Find intercept shift (delta) that matches predicted avg and obs risk
# Root finding used to backsolve
# Expand uniroot brackets until signs differ or ceiling is hit
uniroot_bracket <- function(f, lo = -30, hi = 30, step = 5) {
  flo <- f(lo)
  fhi <- f(hi)
  while (sign(flo) == sign(fhi)) {
    lo <- lo - step
    hi <- hi + step
    flo <- f(lo)
    fhi <- f(hi)
    k <- k + 1
  }
  if (sign(flo) == sign(fhi)) return(NULL)
  uniroot(f, c(lo, hi))$root
}

# Apply functions to find delta for each draw per trial
find_delta <- function(log_b, a, dt, lp, r_obs, N_arm) {
  
  # Get proportion
  y_prop <- r_obs / N_arm
  
  # Compute p from Gompertz CDF
  f <- function(delta) {
    p <- gompertz_prob(b = exp(log_b + delta), a = a, dt = dt, lp = lp)
    mean(p) - y_prop
  }
  
  # Use root-finding for delta
  out <- uniroot_bracket(f, lo = -30, hi = 30, step = 5)
  
  if (is.null(out)) {
    obj <- function(delta) (f(delta))^2
    optimise(obj, interval = c(-40, 40))$minimum
  } else out
  
}

# Setup for backsolving
S <- nrow(b_logb_mat) # Draws
lp_lst <- b_lp # Participant LPs
dt_lst <- b_dt_lst # Time
delta_mat <- matrix(NA_real_, nrow = S, ncol = I) # Matrix to store deltas

# Backsolve to find delta for every draw per trial
for (i in 1:I) {
  lp_i <- lp_lst[[i]]
  dt_i <- dt_lst[[i]]
  r_i <- r[i]
  N_i <- N[i]
  
  # Progress notification
  if (i %% 10 == 1 || i == I) message(sprintf("Arm %d/%d", i, I))
  
  for (s in 1:S) {
    delta_mat[s, i] <- find_delta(
      log_b = b_logb_mat[s, i],
      a = b_a_mat[s, i],
      dt = dt_i,
      lp = lp_i,
      r_obs = r_i,
      N_arm = N_i
    )
  }
}

mean(is.na(delta_mat[,i])) # 0 NAs

# Save delta matrix
saveRDS(delta_mat, "outputs/conduct_model/citl/delta_draws.rds")

# Import to save runtime
delta_mat <- readRDS("outputs/conduct_model/citl/delta_draws.rds")

# Get pre- and post-CITL average risk and predicted counts
avg_risk_pre  <- matrix(NA_real_, nrow = S, ncol = I)
pred_count_pre  <- matrix(NA_real_, nrow = S, ncol = I)

avg_risk_post  <- matrix(NA_real_, nrow = S, ncol = I)
pred_count_post  <- matrix(NA_real_, nrow = S, ncol = I)

for (i in 1:I) {
  lp_i <- b_lp[[i]]
  dt_i <- b_dt_lst[[i]]
  
  # Progress notification
  if (i %% 10 == 1 || i == I) message(sprintf("Arm %d/%d", i, I))
  
  for (s in 1:S) {
    p_pre <- gompertz_prob(
      b = exp(b_logb_mat[s, i]), a = b_a_mat[s, i], 
      dt = dt_i, lp = lp_i
    )
    p_post <- gompertz_prob(
      b = exp(b_logb_mat[s, i] + delta_mat[s, i]), a = b_a_mat[s, i], 
      dt = dt_i, lp = lp_i
    )
    avg_risk_pre[s, i]  <- mean(p_pre)
    avg_risk_post[s, i] <- mean(p_post)
    pred_count_pre[s, i]  <- N[i] * avg_risk_pre[s, i]
    pred_count_post[s, i] <- N[i] * avg_risk_post[s, i]
  }
}

# Summarise per arm
calib_sum <- tibble(
  arm = seq_len(I),
  N = N,
  r_obs = r,
  r_pred = apply(pred_count_post, 2, \(x) median(x, na.rm=TRUE)),
  r_citl = apply(pred_count_post, 2, \(x) median(x, na.rm=TRUE)),
  risk_obs = r/N,
  risk_pred = apply(avg_risk_pre, 2, \(x) median(x, na.rm=TRUE)),
  risk_citl = apply(avg_risk_post, 2, \(x) median(x, na.rm=TRUE)),
  delta_med     = apply(delta_mat, 2, \(x) median(x, na.rm=TRUE)),
  delta_lo      = apply(delta_mat, 2, \(x) quantile(x, 0.025, na.rm=TRUE)),
  delta_hi      = apply(delta_mat, 2, \(x) quantile(x, 0.975, na.rm=TRUE))
)

# Save summary
write_csv(calib_sum, "outputs/conduct_model/citl/delta_risk_sum.csv")

# Check if backsolving hit targets
# At or very close to zero a good sign
res <- avg_risk_post - rep(r / N, nrow(avg_risk_post))
summary(as.numeric(res))

# 64 Sobol
# New model w/o duration, categorical sample size and sites
# Min.     1st Qu.      Median        Mean     3rd Qu.        Max. 
# -0.88091589 -0.08704349 -0.00000009 -0.00000028  0.08704311  0.88092033 

# 64 Sobol
# Min.      1st Qu.   Median   Mean     3rd Qu.  Max. 
# -0.894497 -0.078729 0.009574 0.026206 0.104049 0.985915

# 128 Sobol (very little diff, stick to 64)
# Min.      1st Qu.   Median   Mean     3rd Qu.  Max. 
# -0.894496 -0.078727  0.009573  0.026123  0.104005  0.985915

# Plot
plot_dat <- tibble(
  arm = rep(seq_len(I), 2),
  obs = rep(r / N, 2),
  fit = c(
    apply(avg_risk_pre, 2, median, na.rm = TRUE),
    apply(avg_risk_post, 2, median, na.rm = TRUE)
  ),
  stage = rep(c("pre", "post"), each = I)
)

plot_citl <- plot_dat %>% 
  ggplot(aes(x = obs, y = fit, colour = stage)) +
  geom_abline(linetype = 2) +
  geom_point() +
  labs(
    x = "Observed attrition risk",
    y = "Predicted attrition risk"
  ) +
  theme_classic()

plot_citl

ggsave(
  "outputs/conduct_model/citl/calib_plot.png",
  plot_citl,
  width = 6,
  height = 6,
  units = "in"
)
