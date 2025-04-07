
# Run stan model and inspect outputs
source("scripts/00_functions_packages_etc.R")
a_imp_df <- read_csv("processed_data/df_stan_n340.csv")

# Get stan data
c_stan_data <- list(
  I = nrow(b_fmt2),
  r = b_fmt2$arm_n_attr,
  N = b_fmt2$arm_n,
  dt = b_fmt2$duration,
  arm = as.integer(b_fmt2$type),
  var1 = as.integer(b_fmt2$n_arms),
  var2 = as.vector(scale(b_fmt2$n)),
  var3 = as.vector(scale(b_fmt2$n_facilities))
)

glimpse(c_stan_data)

# Compile stan model
# d_compile <- stan_model("Scripts/03_stan_mdl.stan")
d_compile_1var <- stan_model("Scripts/03_stan_addVars.stan")

# Fit model
d_fit <- sampling(
  d_compile_1var,
  data = c_stan_data,
  warmup = 1000,
  iter = 4000,
  chains = 4,
  seed = 123
)

# Check
traceplot(d_fit, pars = c("a_raw", "beta_arm", "mu1"))
