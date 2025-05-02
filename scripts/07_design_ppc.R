
# Posterior predictive checks for design model
source("scripts/00_functions_packages_etc.R")

# Import data for predictive checks
a_imp_df <- readRDS("processed_data/design_stan_data.rds")

# Import model fit from MARS
a_fit <- readRDS("outputs/design_fit.rds")

# Get summarised draws
a_draws <- as_draws_df(a_fit)
a_sum_draws <- summarise_draws(a_draws)

# Convergence ## ---------------------------------------------------------------

# Split draw into a separate column
b_sum <- a_sum_draws %>% 
  mutate(variable = gsub("\\[[0-9,]+\\]", "", variable)) %>% 
  group_by(variable) %>% 
  mutate(draw = case_when(n() == 1 ~ NA, TRUE ~ row_number())) %>% 
  select(variable, draw, everything()) %>% 
  ungroup()

b_sum_conv <- b_sum %>% 
  group_by(variable) %>% 
  reframe(
    n = n(),
    mean_rhat = mean(rhat, na.rm = TRUE),
    mean_ess = mean(ess_bulk, na.rm = TRUE)
  )

# Trace plots - global parameters (e.g. shrinkage, random effects)
trace_global <- traceplot(
  a_fit, 
  pars = c(
    "lp__", "mean_rate", "sd_rate",
    "sigma_a", "a_mean", "sigma_mu1", "mu1_mean",
    "tau", "lambda", "shrinkage_effective_params"
  )
)

# Trace plots - latent mean for beta shrinkage
trace_z <- traceplot(a_fit, pars = "z")

# Trace plots - local shrinkage
trace_lambda <- traceplot(a_fit, pars = "lambda")

# Trace plots - betas
trace_beta <- traceplot(a_fit, pars = c("beta_arm", "beta"))

# Trace plots - first 50 of primary parameters
trace_a <- traceplot(a_fit, pars = c("beta_arm", "beta"))

# Posterior checks -------------------------------------------------------------

# Compare distribution of attrition counts
plot_ppc <- ppc_dens_overlay(a_imp_df$r, as.matrix(a_fit, pars = "r_new")[1:361, ])
plot_ppc

plot_int <- ppc_intervals(a_imp_df$r, as.matrix(a_fit, pars = "r_new")[1:361, ])
plot_int

# Check pairs
pos_array <- as.array(a_fit)
plot_pairs <- mcmc_pairs(pos_array, pars = c("beta"))
plot_pairs

plot_stat <- ppc_stat(a_imp_df$r, as.matrix(a_fit, pars = "r_new")[1:369, ])
plot_stat

ggsave(
  "outputs/stat.jpg",
  plot_stat,
  width = 6,
  height = 4,
  unit = "in"
)

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
