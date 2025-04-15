
# Inspect outputs
source("scripts/00_functions_packages_etc.R")

# Import model fit from MARS
a_fit <- readRDS("outputs/fit.rds")

# Import counts
a_imp_df <- readRDS("processed_data/stan_dat.rds")

# Check visual convergence
plot_trace <- traceplot(a_fit, pars = c("a_raw", "beta_arm", "mu1"))

ggsave(
  "outputs/trace.png",
  plot_trace,
  width = 8,
  height = 4,
  unit = "in"
)

# Check pairs
pos_array <- as.array(a_fit)
plot_pairs <- mcmc_pairs(pos_array, pars = c("a_raw", "mu1", "beta_arm"))

ggsave(
  "outputs/pairs.png",
  plot_pairs,
  width = 6,
  height = 6,
  unit = "in"
)

# Check counts
plot_ppc <- ppc_dens_overlay(a_imp_df$r, as.matrix(a_fit, pars = "r_new")[1:375, ])

ggsave(
  "outputs/ppc.jpg",
  plot_ppc,
  width = 6,
  height = 4,
  unit = "in"
)

# Check rhat
summary(a_fit)$summary[c("a_raw", "beta_arm", "mu1"),"Rhat"]
summary(a_fit)$summary[c("a_raw", "beta_arm", "mu1"),"n_eff"]

# Check draws
mcmc_hist(as_draws_df(a_fit), pars = "a_raw")
mcmc_hist(as_draws_df(a_fit), pars = "beta_arm")
mcmc_hist(as_draws_df(a_fit), pars = "mu1")
mcmc_hist(as_draws_df(a_fit), pars = "a")
mcmc_hist(as_draws_df(a_fit), pars = "b[1]")
mcmc_hist(as_draws_df(a_fit), pars = "rate")

# Compare observed and posterior
r_new <- extract(a_fit, pars = "r_new")$r_new
r_obs <- a_imp_df$r
ppc_dens_overlay(r_obs, r_new[1:338, ])
ppc_stat(y = r_obs, yrep = r_new[1:338, ])
ppc_intervals(y = r_obs, yrep = r_new[1:338, ])

# Covariate PPCs
covariate_names <- colnames(a_imp_df[["X"]])

# Get plots
ppc_plots <- plot_subgroup_ppcs(
  y = a_imp_df$r,
  yrep = extract(a_fit, pars = "r_new")$r_new,
  covariate_df = as.data.frame(a_imp_df[["X"]]),
  covariate_names = covariate_names
)

# Print plots
for (name in names(ppc_plots)) {
  print(ppc_plots[[name]])
}

# Extract posterior samples
posterior_samples <- as_draws_df(a_fit)
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

# Summarise gompertz parameters
sum_params <- posterior_summary %>% 
  filter(grepl("^a\\[|^b\\[", variable)) %>% 
  mutate(
    variable = gsub("\\[|\\]", "", variable),
    variable = gsub("[0-9]", "", variable)
  ) %>% 
  group_by(variable) %>% 
  reframe(mean = mean(measures.mean))

# Save
saveRDS(posterior_samples, "outputs/posterior_samples.rds")
saveRDS(posterior_summary, "outputs/posterior_sum.rds")
write_csv(posterior_summary, "outputs/posterior_sum.csv")
