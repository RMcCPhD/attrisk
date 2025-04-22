
# Inspect outputs
source("scripts/00_functions_packages_etc.R")

library(shinystan)

# Import model fit from MARS
a_fit <- readRDS("outputs/fit.rds")

# shinystan::launch_shinystan(a_fit)

# Import counts
a_imp_df <- readRDS("processed_data/stan_dat.rds")

# Check pairs
pos_array <- as.array(a_fit)
plot_pairs <- mcmc_pairs(pos_array, pars = c("b"))

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
