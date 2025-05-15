
source("mars/design_model_mk4/scripts/00_config_functions_packages.R")

# Import data, model fit and posterior summaries
a_imp_df <- readRDS("mars/design_model_mk4/stan_dat.rds")
a_imp_mdl <- readRDS("mars/design_model_mk4/outputs/design4_fit.rds")
a_imp_sum <- readRDS("mars/design_model_mk4/outputs/design4_sum_draws.rds")

# LOO invalid because it assumes each is independent
# Not always the case, some trials have two arms (ac and plc)

# Compute WAIC for model
# a_draws <- as.matrix(a_imp_mdl, pars = "log_lik")
# loo::waic(a_draws)

# LOO
# loo_result <- loo(a_draws)
# plot(loo_result)
# which(pareto_k_ids(loo_result) > 0.7)

# Get beta names
b_beta_nm <- colnames(a_imp_df$X)

# Evaluate predictor significance for betas
b_betas <- a_imp_sum %>% 
  filter(grepl("^beta\\[|beta_arm", variable)) %>% 
  mutate(
    beta_name = factor(b_beta_nm, levels = b_beta_nm)
  )

# Save to compare to no shrinkage
# write_csv(b_betas, "mars/design_model_mk4/betas.csv")

plot_b_betas <- b_betas %>% 
  ggplot(aes(x = mean, xmin = `2.5%`, xmax = `97.5%`, y = fct_rev(beta_name))) +
  geom_point() +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, colour = "red") +
  geom_vline(xintercept = 0.05, colour = "orange", linetype = "dashed") +
  geom_vline(xintercept = -0.05, colour = "orange", linetype = "dashed") +
  theme_bw() +
  scale_x_continuous(n.breaks = 20) +
  labs(x = "mean (95% credible intervals) effect", y = "covariate")

ggsave(
  "mars/design_model_mk4/outputs/betas.png",
  plot_b_betas,
  width = 8,
  height = 8,
  units = "in"
)

# Convergence ## ---------------------------------------------------------------

# Import convergence summary
c_imp_sum <- read_csv("mars/design_model_mk4/outputs/sum_rhat_ess.csv")

# Check for divergences
check_hmc_diagnostics(a_imp_mdl) # Nothing abnormal, i.e. 0 divergences

# Trace plots - global parameters (e.g. shrinkage, random effects)
trace_global <- traceplot(
  a_imp_mdl, 
  pars = c(
    "lp__", "mean_rate", "sd_rate",
    "tau_a", "mu_a", "tau_mu1", "mu_mu1",
    "tau_beta"
  )
)

ggsave(
  "mars/design_model_mk4/outputs/convergence/trace_global.png", 
  trace_global, 
  width = 8, 
  height = 8,
  units = "in"
)

# Trace plots - latent mean for beta shrinkage
trace_z <- traceplot(a_imp_mdl, pars = "z_beta")

ggsave(
  "mars/design_model_mk4/outputs/convergence/trace_shrinkage_coef.png", 
  trace_z, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots - local shrinkage
trace_lambda <- traceplot(a_imp_mdl, pars = "lambda_beta")

ggsave(
  "mars/design_model_mk4/outputs/convergence/trace_shrinkage_lambda.png", 
  trace_lambda, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots - betas
trace_beta <- traceplot(a_imp_mdl, pars = "beta")

ggsave(
  "mars/design_model_mk4/outputs/convergence/trace_beta.png", 
  trace_beta, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots - first 50 of primary parameters
trace_a <- traceplot(a_imp_mdl, pars = paste0("a[", 1:56, "]"))
trace_b <- traceplot(a_imp_mdl, pars = paste0("b[", 1:56, "]"))

ggsave(
  "mars/design_model_mk4/outputs/convergence/trace_a.png", 
  trace_a, 
  width = 15, 
  height = 10,
  units = "in"
)

ggsave(
  "mars/design_model_mk4/outputs/convergence/trace_b.png", 
  trace_b, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots - baseline log-hazard rate (mu1)
trace_mu1 <- traceplot(a_imp_mdl, pars = paste0("mu1[", 1:56, "]"))

ggsave(
  "mars/design_model_mk4/outputs/convergence/trace_mu1.png", 
  trace_mu1, 
  width = 15, 
  height = 10,
  units = "in"
)

# Posterior checks -------------------------------------------------------------

# Compare distribution of attrition counts
ppc_dens <- ppc_dens_overlay(a_imp_df$r, as.matrix(a_imp_mdl, pars = "r_new")[1:361, ])
ppc_int <- ppc_intervals(a_imp_df$r, as.matrix(a_imp_mdl, pars = "r_new")[1:361, ])

ggsave(
  "mars/design_model_mk4/outputs/posterior_checks/ppc_density.jpg", 
  ppc_dens, 
  width = 8, 
  height = 8,
  units = "in"
)

ggsave(
  "mars/design_model_mk4/outputs/posterior_checks/ppc_intervals.jpg", 
  ppc_int, 
  width = 8, 
  height = 8,
  units = "in"
)

# Residuals
ppc_ppe_mean <- ppc_stat(
  y = a_imp_df$r, 
  yrep = as.matrix(a_imp_mdl, pars = "r_new")[1:361, ], 
  stat = "mean"
)

ppc_ppe_sd <- ppc_stat(
  y = a_imp_df$r, 
  yrep = as.matrix(a_imp_mdl, pars = "r_new")[1:361, ], 
  stat = "sd"
)

ggsave(
  "mars/design_model_mk4/outputs/posterior_checks/ppc_pred_err_mean.jpg", 
  ppc_ppe_mean, 
  width = 8, 
  height = 8,
  units = "in"
)

ggsave(
  "mars/design_model_mk4/outputs/posterior_checks/ppc_pred_err_sd.jpg", 
  ppc_ppe_sd, 
  width = 8, 
  height = 8,
  units = "in"
)

# Pairwise correlation - coefficients
draws_df <- as_draws_df(a_imp_mdl)
ppc_cor_beta <- draws_df[, grepl("^beta", names(draws_df))]

ppc_cor_beta_df <- as.data.frame(cor(ppc_cor_beta)) %>% 
  rownames_to_column(var = "row") %>% 
  pivot_longer(-row, names_to = "column", values_to = "cor") %>% 
  left_join(b_betas %>% select(row = variable, beta_name)) %>% 
  rename(row_name = beta_name) %>% 
  mutate(column_name = case_when(row == column ~ row_name)) %>% 
  group_by(column) %>% 
  fill(column_name, .direction = "downup") %>% 
  ungroup()

ppc_cor_beta_plot <- ppc_cor_beta_df %>% 
  ggplot(aes(x = row_name, y = column_name, fill = cor)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Beta rows", y = "Beta columns", fill = "R")

ggsave(
  "design_model_mk4/outputs/ppc_cor_beta.jpg",
  ppc_cor_beta_plot,
  width = 8,
  height = 5,
  unit = "in"
)

# Covariate PPCs
covariate_names <- colnames(a_imp_df[["X"]])

# Get plots
ppc_plots <- plot_subgroup_ppcs(
  y = a_imp_df$r,
  yrep = extract(a_imp_mdl, pars = "r_new")$r_new,
  covariate_df = as.data.frame(a_imp_df[["X"]]),
  covariate_names = covariate_names
)

# Print plots
for (name in names(ppc_plots)) {
  print(name)
  save_name <- paste0("ppc_covar_", name, ".jpg")
  path <- paste0("mars/design_model_mk4/outputs/", save_name)
  plot <- ppc_plots[[name]]
  ggsave(path, plot, width = 6, height = 4, unit = "in")
  print(plot)
}
