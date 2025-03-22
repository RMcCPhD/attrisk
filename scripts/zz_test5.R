
source("scripts/00_functions_packages_etc.R")

# Import data - attrition counts, design variables
imp_trials <- read_rds("processed_data/df_plc_ac.rds")
glimpse(imp_trials)

# Take key variables
imp_trials <- imp_trials %>% 
  select(nct_id:duration)

# Create a mock numeric variable
# Remove one trial with 100% attrition (likely a count error)
imp_trials_mock <- imp_trials %>% 
  mutate(mock_age = rnorm(924, mean = 52.5, sd = 10.6)) %>% 
  filter(!is.na(duration))

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

# Summarise parameters
fit_sum <- summarize_draws(fit)

fit_sum %>% 
  filter(grepl("a\\[|b\\[", variable)) %>% 
  mutate(
    variable = gsub("\\[|\\]", "", variable),
    variable = gsub("[0-9]", "", variable)
  ) %>% 
  group_by(variable) %>% 
  reframe(mean = mean(mean))

fit_sum2 <- fit_sum %>% 
  filter(grepl("a\\[|b\\[", variable)) %>% 
  mutate(
    iter = str_extract(variable, "(?<=\\[)\\d+(?=\\])"),
    iter = as.numeric(iter),
    variable = gsub("\\[|\\]", "", variable),
    variable = gsub("[0-9]", "", variable)
  ) %>% 
  arrange(iter, variable)

# write_csv(fit_sum2, "scratch_data/stan_test_outputs.csv")

draw_post <- as_draws_df(fit) %>% 
  pivot_longer(-c(lp__:.draw), names_to = "parameter", values_to = "value")

draw_post %>% select(a, b) %>% cov()

fit_sum2 %>% 
  group_by(iter, variable) %>% 
  reframe(mean = mean(mean)) %>% 
  pivot_wider(names_from = "variable", values_from = "mean")

# Compare observed vs predicted counts
r <- extract(fit, pars = "r_new")$r_new
obs_r <- stan_data$r

ppc_dens_overlay(obs_r, r[1:207, ])
ppc_intervals(y = obs_r, yrep = r[1:207, ])# Check divergent transitions
check_hmc_diagnostics(fit)

# Check trace plots
traceplot(fit, pars = c("a_raw", "beta", "b_mu", "u[1]", "u[2]"))
