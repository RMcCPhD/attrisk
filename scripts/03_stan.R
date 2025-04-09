
# Run stan model and inspect outputs
source("scripts/00_functions_packages_etc.R")

a_imp_df <- readRDS("processed_data/df_stan_ready.rds") %>% rename(b_dbl = b_dlb)
a_imp_df %>% count(trial_id)

# Setup stan dataset (340 trials, 367 rows, 41 vars incl. arm)
c_stan_data <- list(
  I = nrow(a_imp_df),
  r = a_imp_df$arm_n_attr,
  N = a_imp_df$arm_n,
  dt = a_imp_df$duration,
  arm = as.integer(a_imp_df$type),
  var1 = as.integer(a_imp_df$p3),
  var2 = as.integer(a_imp_df$p4),
  var3 = as.integer(a_imp_df$po),
  var4 = as.integer(a_imp_df$b_dbl),
  var5 = as.integer(a_imp_df$b_ol),
  var6 = as.integer(a_imp_df$b_q),
  var7 = as.integer(a_imp_df$b_t),
  var8 = as.integer(a_imp_df$reg_1),
  var9 = as.integer(a_imp_df$reg_2),
  var10 = as.integer(a_imp_df$reg_3),
  var11 = as.integer(a_imp_df$reg_4),
  var12 = as.integer(a_imp_df$reg_5),
  var13 = as.integer(a_imp_df$reg_6),
  var14 = as.integer(a_imp_df$reg_am),
  var15 = as.integer(a_imp_df$reg_eur),
  var16 = as.integer(a_imp_df$reg_pac),
  var17 = as.integer(a_imp_df$reg_afr),
  var18 = as.integer(a_imp_df$reg_sea),
  var19 = as.integer(a_imp_df$reg_med),
  var20 = as.integer(a_imp_df$n_250),
  var21 = as.integer(a_imp_df$n_500),
  var22 = as.integer(a_imp_df$n_750),
  var23 = as.integer(a_imp_df$n_1000),
  var24 = as.integer(a_imp_df$n_abv),
  var25 = as.integer(a_imp_df$arm_2),
  var26 = as.integer(a_imp_df$arm_3),
  var27 = as.integer(a_imp_df$arm_abv),
  var28 = as.integer(a_imp_df$dur_6m),
  var29 = as.integer(a_imp_df$dur_1yr),
  var30 = as.integer(a_imp_df$dur_abv),
  var31 = as.integer(a_imp_df$pom_1),
  var32 = as.integer(a_imp_df$pom_2),
  var33 = as.integer(a_imp_df$pom_3),
  var34 = as.integer(a_imp_df$pom_4),
  var35 = as.integer(a_imp_df$pom_5),
  var36 = as.integer(a_imp_df$pom_6),
  var37 = as.integer(a_imp_df$site_50),
  var38 = as.integer(a_imp_df$site_100),
  var39 = as.integer(a_imp_df$site_150),
  var40 = as.integer(a_imp_df$site_abv)
)

glimpse(c_stan_data)

# Compile stan model
d_compile <- stan_model("Scripts/03_model.stan")

# Fit model
d_fit <- sampling(
  d_compile,
  data = c_stan_data,
  warmup = 1000,
  iter = 3000,
  chains = 4,
  seed = 123
)

# Check visual convergence
traceplot(d_fit, pars = c("a_raw", "beta_arm", "mu1"))

# Check rhat
summary(d_fit)$summary[c("a_raw", "beta_arm", "mu1"),"Rhat"]
summary(d_fit)$summary[c("a_raw", "beta_arm", "mu1"),"n_eff"]

# Check draws
mcmc_hist(as_draws_df(d_fit), pars = "a_raw")
mcmc_hist(as_draws_df(d_fit), pars = "beta_arm")
mcmc_hist(as_draws_df(d_fit), pars = "mu1")
mcmc_hist(as_draws_df(d_fit), pars = "a")
mcmc_hist(as_draws_df(d_fit), pars = "b[1]")
mcmc_hist(as_draws_df(d_fit), pars = "rate")

# Compare observed and posterior
r_new <- extract(d_fit, pars = "r_new")$r_new
r_obs <- c_stan_data$r
ppc_dens_overlay(r_obs, r_new[1:338, ])
ppc_stat(y = r_obs, yrep = r_new[1:338, ])
ppc_intervals(y = r_obs, yrep = r_new[1:338, ])

# Covariate PPCs
covariate_names <- names(select(a_imp_df, -c(trial_id:duration)))

# Get plots
ppc_plots <- plot_subgroup_ppcs(
  y = a_imp_df$arm_n_attr,
  yrep = extract(d_fit, pars = "r_new")$r_new,
  covariate_df = a_imp_df,
  covariate_names = covariate_names
)

# Print plots
for (name in names(ppc_plots)) {
  print(ppc_plots[[name]])
}

# Extract posterior samples
posterior_samples <- as_draws_df(d_fit)
posterior_summary <- summarise_draws(posterior_samples, default_summary_measures())

posterior_summary <- summarise_draws(
  posterior_samples,
  measures = list(
    mean = mean,
    median = median,
    q5 = ~quantile(.x, 0.05),
    q95 = ~quantile(.x, 0.95)
  )
)

# Save
saveRDS(posterior_samples, "outputs/posterior_samples.rds")
saveRDS(posterior_summary, "outputs/posterior_sum.rds")
write_csv(posterior_summary, "outputs/posterior_sum.csv")
