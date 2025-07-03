
source("mars/design_model_mk5/scripts/00_config_functions_packages.R")

a_imp_df <- readRDS("mars/design_model_mk5/stan_dat.rds")

# Add variables needed for aggregate logistic formula
b_df <- tibble(
  trial_id = a_imp_full$trial_id,
  arm_n_attr = a_imp_df$r,
  arm_n_compl = a_imp_df$N - a_imp_df$r,
  dt = a_imp_df$dt,
  log_dt = log(a_imp_df$dt)
)

# Test with intercept-only model
b_formula <- as.formula(paste("arm_n_attr | trials(arm_n_attr + arm_n_compl) ~ 1 + (1 | trial_id) + offset(log_dt)"))

b_fit <- brm(
  formula = b_formula,
  data = b_df,
  family = binomial(),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE)
)

summary(b_fit)
b_draws <- as_draws_df(b_fit)

# Posterior checks -------------------------------------------------------------

# Cross-validation - looks good, all k < 0.7
loo(b_fit, moment_match = TRUE)
bayes_R2(b_fit) # 0.989

# Trace plot and recovery
trace <- mcmc_trace(c_draws, pars = c("b_Intercept", "Intercept", "lprior", "lp__"))
r_overall <- pp_check(b_fit, type = "dens_overlay", ndraws = 361)
r_ind <- pp_check(b_fit, type = "intervals")

ggsave(
  "C://Users/ryanm/Desktop/meeting/brms_intercept_only/trace.jpg",
  trace,
  height = 8,
  width = 8,
  units = "in"
)

ggsave(
  "C://Users/ryanm/Desktop/meeting/brms_intercept_only/r_recovery_overall.jpg",
  r_overall,
  height = 8,
  width = 8,
  units = "in"
)

ggsave(
  "C://Users/ryanm/Desktop/meeting/brms_intercept_only/r_recovery_individual.jpg",
  r_ind,
  height = 8,
  width = 8,
  units = "in"
)

sum_draws <- summarise_draws(b_draws)
write_csv(sum_draws, "C://Users/ryanm/Desktop/meeting/brms_intercept_only/sum_draws.csv")
