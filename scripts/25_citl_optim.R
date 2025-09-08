
# Solve intercept shift using 1-D backsolving
# Conduct model becomes two-stage
# 1. Fit pre-CITL model to obtain posterior draws (uncalibrated rate)
# 2. Backsolve to adjust intercept to match arm-level average risk

source("scripts/00_config.R")
library(tidyverse)
library(posterior)
library(matrixStats)

# Import CITL data and flattened posterior arrays
a_imp_dat <- readRDS("mars/citl_model/latest_stan_citl.rds")
a_imp_mu1 <- readRDS("mars/citl_model/citl_mu1.rds")
a_imp_beta <- readRDS("mars/citl_model/citl_beta.rds")
a_imp_error <- readRDS("mars/citl_model/citl_error.rds")

# Get dims
I <- a_imp_dat$I
K <- a_imp_dat$K
J <- a_imp_dat$J
total_draws <- nrow(a_imp_mu1)

# Configure sampling
set.seed(123)
D_target <- 4000L
out_dir <- sprintf("mars/citl_model/citl_backsolve_%03d", D_target)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Pick draws evenly across total
ids <- unique(as.integer(round(seq(1, total_draws, length.out = D_target))))

# Compute arm-level HRs
hr_all <- as.vector(exp(a_imp_dat$Z %*% a_imp_dat$pB))
hr_lst <- split(hr_all, rep(seq_len(I), each = J))

# Pull out parameters
r <- a_imp_dat$r
N <- a_imp_dat$N
dt <- a_imp_dat$dt
lp <- a_imp_dat$lp
beta_rate <- a_imp_dat$beta_rate
X <- a_imp_dat$X

# Compute average risk for each arm at intercept shift (gamma) for a given
# log-rate and sampled draws
avg_risk <- function(gamma, i, log_b_pre_i, error_i, dt_i, hr_i) {
  log_b_cal <- log_b_pre_i + gamma
  b <- exp(log_b_cal)
  a <- log_b_cal * beta_rate + lp[i] + error_i
  
  # Compute Gompertz constant
  C_i <- if (abs(a) < 1e-4) b * dt_i else (b / a) * expm1(a * dt_i)
  
  # Numerically stable average risk
  mlog <- logSumExp(-hr_i * C_i) - log(J)
  pmin(pmax(1 - exp(mlog), 1e-12), 1 - 1e-12) # Clamp to prevent exact 0/1
}

# Solve gamma for each arm, matching average and observed risk
solve_gamma <- function(i, log_b_pre_i, error_i) {
  target <- pmin(pmax(r[i] / N[i], 1e-12), 1 - 1e-12) # Clamped
  f <- function(g) avg_risk(g, i, log_b_pre_i, error_i, dt[i], hr_lst[[i]]) - target
  
  g0 <- 0
  f0 <- f(g0)
  if (!is.finite(f0) || abs(f0) < 1e-12) return(g0)
  
  # Monotonic increasing bracket root, expand to 60
  step <- 1
  max_expand <- 60
  
  if (f0 > 0) {
    lo <- g0
    for (k in 1:max_expand) {
      lo <- lo - step
      if (f(lo) <= 0) break
      step <- step * 2
    }
    hi <- lo + step
  } else {
    hi <- g0
    for (k in 1:max_expand) {
      hi <- hi + step
      if (f(hi) >= 0) break
      step <- step * 2
    }
    lo <- hi - step
  }
  
  f_lo <- f(lo)
  f_hi <- f(hi)
  if (f_lo * f_hi > 0) return(if (abs(f_lo) < abs(f_hi)) lo else hi)
  uniroot(f, lower = lo, upper = hi, tol = 1e-8)$root
}

# Backsolve per draw
res_index <- tibble(draw_id = ids)
sums <- vector("list", length(ids))

for (j in seq_along(ids)) {
  
  # Setup draws
  d <- ids[j]
  log_b_pre <- as.vector(a_imp_mu1[d, ] + X %*% a_imp_beta[d, ])
  error_d <- a_imp_error[d, ]
  
  # Solve gamma per arm
  gamma_solv <- vapply(
    seq_len(I), 
    function(i) solve_gamma(i, log_b_pre[i], error_d[i]), numeric(1)
  )
  
  # Derived quantities
  delta <- gamma_solv
  haz_mult <- exp(delta)
  
  # Compute risk at solved gamma
  p_solv <- vapply(seq_len(I), function(i) {
    avg_risk(gamma_solv[i], i, log_b_pre[i], error_d[i], dt[i], hr_lst[[i]])
  }, numeric(1))
  
  # Save per-draw outputs
  out <- list(
    draw_id = d,
    gamma = gamma_solv,
    delta = delta,
    haz_mult = haz_mult,
    p_solv = p_solv
  )
  
  saveRDS(out, file = file.path(out_dir, sprintf("citl_backsolve_d%05d.rds", d)))
  
  # Compact summary
  sums[[j]] <- tibble(
    draw_id = d,
    gamma_mean = mean(gamma_solv),
    gamma_sd = sd(gamma_solv),
    p_solv_mean = mean(p_solv),
    p_solve_rmse = sqrt(mean((p_solv - r / N) ^ 2))
  )
  
  if (j %% 10 == 0 || j == length(ids)) {
    message(sprintf("[%d/%d] done draw %d", j, length(ids), d))
  }
}

# Aggregate over draws
files <- list.files(
  out_dir, pattern = "^citl_backsolve_d\\d+\\.rds$", full.names = TRUE
)

# Per arm summary of gamma and hazard multiplier
arm_sum <- map_dfr(files, function(f) {
  obj <- readRDS(f)
  tibble(
    draw_id = obj$draw_id, 
    arm = seq_len(I),
    gamma = obj$gamma, 
    haz_mult = obj$haz_mult
  )
})

gamma_sum <- arm_sum %>% 
  group_by(arm) %>% 
  reframe(
    gamma_mean = mean(gamma),
    gamma_l95  = quantile(gamma, 0.025),
    gamma_u95  = quantile(gamma, 0.975),
    haz_mean   = mean(haz_mult),
    haz_l95    = quantile(haz_mult, 0.025),
    haz_u95    = quantile(haz_mult, 0.975)
  )

# Save backsolving outputs
write_csv(bind_rows(sums), file.path(out_dir, "draw_level_summaries.csv"))
write_csv(gamma_sum, file.path(out_dir, "arm_level_gamma_summary.csv"))

message("Back-solving complete. Outputs in: ", out_dir)

# Sanity check -----------------------------------------------------------------

# 0) Basic structure checks ----------------------------------------------------
stopifnot(length(hr_all) == I * J)
stopifnot(all(vapply(hr_lst, length, integer(1)) == J))
stopifnot(nrow(a_imp_mu1) == nrow(a_imp_beta), nrow(a_imp_mu1) == nrow(a_imp_error))
stopifnot(ncol(a_imp_mu1) == I, ncol(a_imp_beta) == K, ncol(a_imp_error) == I)

# 1) Choose an arm i and stage-1 draw d ---------------------------------------
set.seed(1)
i <- 1                     # pick an arm to inspect (change as needed)
d <- 1                     # pick a stage-1 draw index (change as needed)

# You can also pick a random draw:
# d <- sample.int(nrow(a_imp_mu1), 1)

# 2) Prepare draw-specific pieces for that arm --------------------------------
log_b_pre_vec <- as.vector(a_imp_mu1[d, ] + X %*% a_imp_beta[d, ])
error_vec     <- a_imp_error[d, ]
log_b_pre_i   <- log_b_pre_vec[i]
error_i       <- error_vec[i]
hr_i          <- hr_lst[[i]]
target        <- pmin(pmax(r[i] / N[i], 1e-12), 1 - 1e-12)

# 3) Define the arm-level average-risk function (PROBABILITY, not log) --------
avg_risk_arm_prob <- function(gamma) {
  log_b_cal <- log_b_pre_i + gamma
  b <- exp(log_b_cal)
  a <- beta_rate * log_b_cal + lp[i] + error_i
  # Gompertz constant
  C <- if (abs(a) < 1e-4) b * dt[i] else (b / a) * expm1(a * dt[i])
  # log-mean-exp for stability
  mlog <- matrixStats::logSumExp(-hr_i * C) - log(J)
  # convert to probability and clamp
  pmin(pmax(1 - exp(mlog), 1e-12), 1 - 1e-12)
}

# 4) Visual check: p(g) should be monotone increasing and cross target --------
grid <- seq(-20, 20, by = 2)
pg   <- tibble(g = grid, p = vapply(grid, avg_risk_arm_prob, 0.0))

print(pg, n = Inf)
cat("target r/N:", target, "\n")
cat("p(g) range on grid:", range(pg$p), "\n")
cat("Monotone increasing? ",
    all(diff(pg$p) > -1e-10),  # allow tiny numerical wiggle
    "\n\n")

# Optional quick plot
# ggplot(pg, aes(g, p)) + geom_line() + geom_hline(yintercept = target, linetype = 2)

# 5) Robust bracketing + root solve for gamma* --------------------------------
bracket_mono <- function(f, g0 = 0, step = 1, max_expand = 60) {
  f0 <- f(g0)
  if (!is.finite(f0) || abs(f0) < 1e-12) return(c(g0, g0))
  if (f0 > 0) {  # need smaller gamma
    lo <- g0
    for (k in 1:max_expand) { lo <- lo - step; if (f(lo) <= 0) return(c(lo, g0)); step <- step * 2 }
    return(c(lo, g0))
  } else {       # need larger gamma
    hi <- g0
    for (k in 1:max_expand) { hi <- hi + step; if (f(hi) >= 0) return(c(g0, hi)); step <- step * 2 }
    return(c(g0, hi))
  }
}

f <- function(g) avg_risk_arm_prob(g) - target
br <- bracket_mono(f, g0 = 0)

if (br[1] == br[2]) {
  gamma_hat <- br[1]
  cat("Exact hit at g0; gamma* =", gamma_hat, "\n")
} else {
  f_lo <- f(br[1]); f_hi <- f(br[2])
  if (f_lo * f_hi > 0) {
    gamma_hat <- if (abs(f_lo) < abs(f_hi)) br[1] else br[2]
    cat("Could not straddle root; using boundary gamma* =", gamma_hat, "\n")
  } else {
    gamma_hat <- uniroot(f, lower = br[1], upper = br[2], tol = 1e-8)$root
  }
}

# 6) Validate solution ---------------------------------------------------------
p_hat <- avg_risk_arm_prob(gamma_hat)
cat(sprintf("Solved gamma* = %.6f  =>  p(gamma*) = %.6f  (target = %.6f)  diff = %.3e\n",
            gamma_hat, p_hat, target, p_hat - target))

# 7) If p(g) never reaches target on grid, suggest expanding range ------------
if (target < min(pg$p) || target > max(pg$p)) {
  cat("NOTE: target lies outside p(g) over [-20,20]. Try a wider grid.\n")
}

# 8) (Optional) quick checks for Z stacking & units ----------------------------
# Check that the first arm's hr_i matches slicing hr_all
hr_i_slice <- hr_all[((i - 1) * J + 1):(i * J)]
stopifnot(all.equal(hr_i, hr_i_slice, tolerance = 0, check.attributes = FALSE))

# If p(g) is ~1 for all g or ~0 for all g, inspect dt units and magnitudes of hr_i and beta_rate.


