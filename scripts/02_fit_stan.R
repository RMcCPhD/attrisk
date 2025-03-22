
# Get stan data
stan_data <- list(
  I = nrow(imp_trials_mock),
  r = imp_trials_mock$arm_n_attr,
  N = imp_trials_mock$arm_n,
  dt = imp_trials_mock$duration,
  var1 = imp_trials_mock$mock_age,
  arm = as.integer(imp_trials_mock$flag)
)

str(stan_data)

# Compile model
compile <- stan_model("Scripts/zz_test6.stan")

# Fit model
fit <- sampling(
  compile,
  data = stan_data,
  warmup = 1000,
  iter = 4000,
  chains = 3,
  seed = 123
)

# Check convergence
traceplot(
  fit, 
  pars = c(
    "a_raw", "beta", "b_mu", "u[1]", "u[2]", 
    "a[1]", "a[2]", 
    "b[1]", "b[2]", 
    "rate"
  )
)

check_hmc_diagnostics(fit)

rhat(fit_draws)
ess_bulk(draws)
ess_tail(draws)

# Compare observed vs predicted counts
r <- extract(fit, pars = "r_new")$r_new
obs_r <- stan_data$r
ppc_dens_overlay(obs_r, r[1:207, ])
ppc_intervals(y = obs_r, yrep = r[1:207, ])

# Extract posterior draws
fit_draws <- as_draws_df(fit)
summarise_draws(fit_draws)

draw_vars <- fit_draws %>% select(a_raw, beta, gamma, sigma_u, b_mu)

attrition_probs <- with(draw_vars, {
  a <- a_raw + beta * 45.17855  + gamma * 1
  b <- exp(b_mu + beta * 45.17855 + gamma * 1 + 0.411067961)
  1 - exp(-b/a * (exp(a * 26) - 1))
})

# Estimating hazard function
gomp_draw <- fit_draws %>% 
  select(contains(c("a[", "b["))) %>% 
  pivot_longer(
    everything(),
    names_to = "parameter",
    values_to = "draw"
  ) %>% 
  mutate(
    iter = str_extract(parameter, "(?<=\\[)\\d+(?=\\])"),
    iter = as.integer(iter),
    parameter = if_else(grepl("a", parameter), "a", "b")
  ) %>% 
  pivot_wider(
    names_from = "parameter",
    values_from = "draw"
  )
