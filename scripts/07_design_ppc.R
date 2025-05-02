
# Posterior predictive checks for design model
source("scripts/00_functions_packages_etc.R")

# Import data for predictive checks
a_imp_df <- readRDS("processed_data/trialsub_stan_data.rds")

# Import model fit from MARS
a_fit <- readRDS("outputs/design_model/design_fit.rds")

# Convergence ## ---------------------------------------------------------------

# Summarise rhat and ess
b_conv_sum <- summary(a_fit)$summary[, "Rhat"] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "parameter") %>% 
  rename(rhat = ".") %>% 
  mutate(
    og_parameter = parameter,
    parameter = gsub("\\[|\\]", "", parameter),
    parameter = gsub("[0-9]", "", parameter)
  ) %>% 
  group_by(parameter) %>% 
  reframe(
    n = n(),
    # min_rhat = min(rhat, na.rm = T),
    # max_rhat = max(rhat, na.rm = T),
    mean_rhat = mean(rhat, na.rm = T),
    # med_rhat = median(rhat, na.rm = T)
  ) %>% 
  left_join(
    summary(a_fit)$summary[, "n_eff"] %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "parameter") %>% 
      rename(ess = ".") %>% 
      mutate(
        og_parameter = parameter,
        parameter = gsub("\\[|\\]", "", parameter),
        parameter = gsub("[0-9]", "", parameter)
      ) %>% 
      group_by(parameter) %>% 
      reframe(
        n = n(),
        # min_ess = min(ess, na.rm = T),
        # max_ess = max(ess, na.rm = T),
        mean_ess = mean(ess, na.rm = T),
        # med_ess = median(ess, na.rm = T)
      )
  )

# Trace plots for global parameters
b_conv_sum$parameter[b_conv_sum$n == 1]
b_vars <- c(
  "a_mean", "beta_arm", "lp__", "mean_rate", "mu1_mean", "sd_rate",
  "shrinkage_effective_params", "sigma_a", "sigma_mu1", "tau"
)
b_trace <- mcmc_trace(as.array(a_fit), pars = b_vars)


# Extract posterior samples
posterior_samples <- as_draws_df(a_fit)
posterior_summary <- summarise_draws(posterior_samples, default_summary_measures())

posterior_summary <- summarise_draws(
  posterior_samples,
  measures = list(
    mean = mean,
    median = median,
    q5 = ~quantile(.x, 0.05),
    q97.5 = ~quantile(.x, 0.975)
  )
)

rhat_summary <- summary(a_fit)$summary[, "Rhat"] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "parameter") %>% 
  rename(rhat = ".") %>% 
  mutate(
    parameter = gsub("\\[|\\]", "", parameter),
    parameter = gsub("[0-9]", "", parameter)
  ) %>% 
  group_by(parameter) %>% 
  reframe(
    n = n(),
    min_rhat = min(rhat, na.rm = T),
    max_rhat = max(rhat, na.rm = T),
    mean_rhat = mean(rhat, na.rm = T),
    med_rhat = median(rhat, na.rm = T)
  )

ess_summary <- summary(a_fit)$summary[, "n_eff"] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "parameter") %>% 
  rename(rhat = ".") %>% 
  mutate(
    parameter = gsub("\\[|\\]", "", parameter),
    parameter = gsub("[0-9]", "", parameter)
  ) %>% 
  group_by(parameter) %>% 
  reframe(
    n = n(),
    min_rhat = min(rhat, na.rm = T),
    max_rhat = max(rhat, na.rm = T),
    mean_rhat = mean(rhat, na.rm = T),
    med_rhat = median(rhat, na.rm = T)
  )

# Check pairs
pos_array <- as.array(a_fit)
plot_pairs <- mcmc_pairs(pos_array, pars = c("beta"))

ggsave(
  "outputs/pairs_b.png",
  plot_pairs,
  width = 10,
  height = 10,
  unit = "in"
)

# Check counts
plot_ppc <- ppc_dens_overlay(a_imp_df$r, as.matrix(a_fit, pars = "r_new")[1:369, ])
plot_ppc

ggsave(
  "outputs/ppc.jpg",
  plot_ppc,
  width = 6,
  height = 4,
  unit = "in"
)

plot_stat <- ppc_stat(a_imp_df$r, as.matrix(a_fit, pars = "r_new")[1:369, ])
plot_stat

ggsave(
  "outputs/stat.jpg",
  plot_stat,
  width = 6,
  height = 4,
  unit = "in"
)

plot_int <- ppc_intervals(a_imp_df$r, as.matrix(a_fit, pars = "r_new")[1:369, ])
plot_int

ggsave(
  "outputs/int.jpg",
  plot_int,
  width = 6,
  height = 4,
  unit = "in"
)

# Check draws
mcmc_hist(as_draws_df(a_fit), pars = "mu1")
mcmc_hist(as_draws_df(a_fit), pars = "rate")
mcmc_hist(as_draws_df(a_fit), pars = "a[1]")
mcmc_hist(as_draws_df(a_fit), pars = "beta_arm")
mcmc_hist(as_draws_df(a_fit), pars = "b[1]")

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
  print(name)
  save_name <- paste0("ppc_covar_", name, ".jpg")
  path <- paste0("outputs/", save_name)
  plot <- ppc_plots[[name]]
  ggsave(path, plot, width = 6, height = 4, unit = "in")
  print(plot)
}
