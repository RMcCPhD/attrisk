
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## SUMMARY OF SCRIPT

#' ## 1. TESTING STAN MODEL

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

source("scripts/00_functions_packages_etc.R")

## 1. TESTING STAN MODEL #######################################################

# Import data - attrition counts, design variables
imp_trials <- read_rds("processed_data/df_oattr.rds")
glimpse(imp_trials)

# Remove NA regions for now
imp_trials <- na.omit(imp_trials)

# Get stan data
stan_data <- list(
  N = nrow(imp_trials),
  r = imp_trials$attr,
  n = imp_trials$n,
  dt = imp_trials$duration,
  var1 = imp_trials$number_of_arms,
  var2 = imp_trials$n_regions
)

str(stan_data)

# Compile model
compile <- stan_model("Scripts/zz_test_stan_spec.stan")

# Fit model
fit <- sampling(
  compile,
  data = stan_data,
  iter = 8000,
  chains = 4,
  seed = 123
)

# Summarise
print(
  fit,
  pars = c("mu1", "beta", "a")
)

# Trace plot
stan_trace(fit, pars = c("mu1", "beta", "a"))

# Compare observed vs predicted counts
r <- extract(fit, pars = "r_new")$r_new
obs_r <- stan_data$r
ppc_dens_overlay(obs_r, r[1:100, ])

# Compare test statistics
ppc_stat(obs_r, r, stat = "mean")
ppc_stat(obs_r, r, stat = "sd")

## 1. TESTING STAN MODEL 2 ######################################################

# Import data - attrition counts, design variables
imp_trials <- read_rds("processed_data/df_oattr.rds")
glimpse(imp_trials)

# Take key variables
imp_trials <- imp_trials %>% select(nct_id:duration)

# Create a mock numeric variable
imp_trials_mock <- imp_trials %>% 
  mutate(mock_age = rnorm(207, mean = 52.5, sd = 10.6))

# Get stan data
stan_data <- list(
  I = nrow(imp_trials_mock),
  r = imp_trials_mock$attr,
  N = imp_trials_mock$n,
  dt = imp_trials_mock$duration,
  var1 = imp_trials_mock$mock_age
)

str(stan_data)

# Compile model
compile <- stan_model("Scripts/zz_test3.stan")

# Fit model
fit <- sampling(
  compile,
  data = stan_data,
  iter = 1000,
  chains = 4,
  seed = 123
)

# Summarise
print(
  fit,
  pars = c("mu1", "beta", "a")
)

# Trace plot
stan_trace(fit, pars = c("mu1", "beta", "a"))

# Compare observed vs predicted counts
r <- extract(fit, pars = "r_new")$r_new
obs_r <- stan_data$r
ppc_dens_overlay(obs_r, r[1:100, ])

# Compare test statistics
ppc_stat(obs_r, r, stat = "mean")
ppc_stat(obs_r, r, stat = "sd")

