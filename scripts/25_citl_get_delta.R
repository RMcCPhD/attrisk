
# Find delta (i.e. intercept shift) for CITL
# Pull in posterior draws from conduct model, Sobols and pooled HRs
# Backsolve with uniroot to find delta that aligns observed risk in trials
# with average risk among participants

source("scripts/00_config.R")
library(tidyverse)
library(multinma)
library(posterior)

# Import pre-CITL model, get draws
a_imp_mdl <- readRDS("mars/conduct_model/outputs/conduct_fit.rds")
a_draws <- as_draws(a_imp_mdl)

# Import pre-CITL dataset
a_imp_standf <- readRDS("mars/conduct_model/latest_stan_conduct.rds")

# Import sobol outputs
a_imp_sobol <- readRDS("processed_data/sobol_lst.rds")

# Get posterior draws of log-rate, rate and shape
b_logb <- extract_variable_array(a_draws, variable = "log_b")
b_b <- extract_variable_array(a_draws, variable = "b")
b_a <- extract_variable_array(a_draws, variable = "a")

# Construct DxI matrices from parameter arrays
arr_to_mat <- function(arr) {
  iter <- dim(arr)[1]
  chain <- dim(arr)[2]
  trial <- dim(arr)[3]
  matrix(arr, nrow = iter * chain, ncol = trial)
}

b_logb_mat <- arr_to_mat(b_logb)
b_b_mat <- arr_to_mat(b_b)
b_a_mat <- arr_to_mat(b_a)

# Pull out trials, counts (i.e. target), arm size and time
I <- a_imp_standf$I
r <- a_imp_standf$r
N <- a_imp_standf$N
dt <- a_imp_standf$dt

# Build arm-level participant linear predictors, split by arm
b_sobol_lp <- as.vector(a_imp_sobol$Z %*% a_imp_sobol$loghr)
b_lp <- split(b_sobol_lp, rep(seq_len(I), each = 64L))

# Build time list
b_dt_lst <- lapply(seq_len(I), \(i) rep(dt[i], length(b_lp[[i]])))

# Helper functions
# Compute probability using Gompertz CDF, limit to exp if shape close to 0
gompertz_prob <- function(b, a, dt, lp) {
  if (abs(a) < 1e-4) {
    1 - exp(-b * exp(lp) * dt)
  } else {
    1 - exp(-(b / a) * (exp(a * dt) - 1) * exp(lp))
  }
}

# Expand uniroot brackets until signs differ or ceiling is hit
uniroot_bracket <- function(f, lo = -10, hi = 10, step = 5, max_expand = 6) {
  flo <- f(lo)
  fhi <- f(hi)
  k <- 0
  while (sign(flo) == sign(fhi) && k < max_expand) {
    lo <- lo - step
    hi <- hi + step
    flo <- f(lo)
    fhi <- f(hi)
    k <- k + 1
  }
  if (sign(flo) == sign(fhi)) return(NULL)
  uniroot(f, c(lo, hi))$root
}

# Find delta for d draws
find_delta <- function(log_b, a, dt, lp, r_obs, N_arm) {
  
  # Get proportion
  y_prop <- r_obs / N_arm
  
  # Compute p from Gompertz CDF
  f <- function(delta) {
    p <- gompertz_prob(b = exp(log_b + delta), a = a, dt = dt, lp = lp)
    mean(p) - y_prop
  }
  
  # Check boundary saturation
  f_lo <- f(-30)
  f_hi <- f(30)
  if (is.na(f_lo) || is.na(f_hi)) return(NA_real_)
  if (f_lo > 0 && f_hi > 0) return(NA_real_)
  if (f_lo < 0 && f_hi < 0) return(NA_real_)
  
  # Backsolve to find delta
  out <- uniroot_bracket(f, lo = -10, hi = 10, step = 6, max_expand = 6)
  if (is.null(out)) {
    obj <- function(delta) (f(delta))^2
    optimise(obj, interval = c(-40, 40))$minimum
  } else out
  
}

# Backsolve to find delta for every draw and arm
S <- nrow(b_logb_mat) # Draws
lp_lst <- b_lp # Participant LPs
dt_lst <- b_dt_lst # Time
delta_mat <- matrix(NA_real_, nrow = S, ncol = I) # Matrix to store deltas

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

