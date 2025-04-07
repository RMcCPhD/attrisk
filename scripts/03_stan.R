
# Run stan model and inspect outputs
source("scripts/00_functions_packages_etc.R")
a_imp_df <- read_csv("processed_data/df_stan_n340.csv")

# Categorise variables - 338 removing notApp phase trials
a_cats <- a_imp_df %>% 
  mutate(
    cat_sites = case_when(
      n_facilities <= 50 ~ 0,
      n_facilities >= 50 ~ 1
    ),
    cat_regions = case_when(
      n_regions <= 2 ~ 0,
      TRUE ~ 1
    ),
    cat_duration = case_when(
      duration <= 26 ~ 0,
      duration > 26 ~ 1
    ),
    cat_poms = case_when(
      n_poms == 1 ~ 0,
      n_poms > 1 ~ 1
    )
  ) %>% 
  filter(phase_notapp != 1) %>% 
  select(-phase_notapp)

# Get stan data
c_stan_data <- list(
  I = nrow(a_cats),
  r = a_cats$arm_n_attr,
  N = a_cats$arm_n,
  dt = a_cats$duration,
  arm = as.integer(a_cats$type),
  var1 = as.integer(a_cats$n_arms),
  var2 = as.vector(scale(a_cats$n)),
  var3 = as.integer(a_cats$cat_regions),
  var4 = as.integer(a_cats$cat_sites),
  var5 = as.integer(a_cats$phase_3),
  var6 = as.integer(a_cats$phase_4),
  var7 = as.integer(a_cats$region_of_the_americas),
  var8 = as.integer(a_cats$european_region),
  var9 = as.integer(a_cats$western_pacific_region),
  var10 = as.integer(a_cats$african_region),
  var11 = as.integer(a_cats$south_east_asia_region),
  var12 = as.integer(a_cats$eastern_mediterranean_region),
  var13 = as.integer(a_cats$double),
  var14 = as.integer(a_cats$open_label),
  var15 = as.integer(a_cats$quadruple),
  var16 = as.integer(a_cats$triple),
  var17 = as.integer(a_cats$cat_duration),
  var18 = as.integer(a_cats$cat_poms)
)

glimpse(c_stan_data)

# Compile stan model
# d_compile_1var <- stan_model("Scripts/03_stan_addVars.stan")
d_compile_reg <- stan_model("Scripts/03_stan_mdl.stan")

# Fit model
d_fit <- sampling(
  d_compile_reg,
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
covariate_names <- c(
  "type", "n_arms", "n", "cat_regions", "cat_sites", "phase_3", "phase_4",
  "region_of_the_americas", "european_region", "western_pacific_region", 
  "african_region", "south_east_asia_region", "eastern_mediterranean_region", 
  "double", "open_label", "quadruple", "triple", "cat_duration", 
  "cat_poms"
)

# Get plots
ppc_plots <- plot_subgroup_ppcs(
  y = a_cats$arm_n_attr,
  yrep = extract(d_fit, pars = "r_new")$r_new,
  covariate_df = a_cats,
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
