
# Inspect outputs
source("scripts/00_functions_packages_etc.R")

# Import model fit from MARS
a_fit <- readRDS("outputs/fit2.rds")

# Import counts
a_imp_df <- readRDS("processed_data/stan_dat.rds")

# Summarise
b_sum <- summary(a_fit)$summary %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "parameter")

# saveRDS(b_sum, "outputs/fit2_sum.rds")

# Check visual convergence
plot_trace <- traceplot(a_fit, pars = c("sigma_a", "beta_arm", "mu1", "rate"))
plot_trace

ggsave(
  "outputs/trace.png",
  plot_trace,
  width = 8,
  height = 4,
  unit = "in"
)

plot_trace2 <- traceplot(a_fit, pars = c("beta"))
plot_trace2

ggsave(
  "outputs/trace_betas.png",
  plot_trace2,
  width = 8,
  height = 4,
  unit = "in"
)

# Check Rhat
summary(b_sum[grepl("beta", b_sum$parameter),]$Rhat)
summary(b_sum[grepl("a", b_sum$parameter),]$Rhat)
summary(b_sum[grepl("b", b_sum$parameter),]$Rhat)

summary(a_fit)$summary[c("sigma_a", "beta_arm", "mu1", "rate"), "Rhat"]
summary(a_fit)$summary[c("sigma_a", "beta_arm", "mu1", "rate"), "n_eff"]

# Check pairs
# pos_array <- as.array(a_fit)
# plot_pairs <- mcmc_pairs(pos_array, pars = c("mu1", "beta_arm", "rate"))
# plot_pairs
# 
# ggsave(
#   "outputs/pairs.png",
#   plot_pairs,
#   width = 6,
#   height = 6,
#   unit = "in"
# )

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
write_csv(sum_params, "outputs/gomp_sum.csv")
