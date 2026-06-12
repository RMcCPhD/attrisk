
# Posterior checks for trial-level sub-model
source("scripts/00_functions_packages_etc.R")

# Import data for predictive checks
a_imp_df <- readRDS("processed_data/trialsub_stan_data.rds")

# Import model fit from MARS
a_fit <- readRDS("outputs/trial_submodel_fit.rds")

# Get summarised draws
a_draws <- as_draws_df(a_fit)
a_sum_draws <- summarise_draws(a_draws)

# Get array for ppc
a_array <- as.array(a_fit)

# Convergence ## ---------------------------------------------------------------

# Check for divergences
check_hmc_diagnostics(a_fit) # Nothing abnormal, i.e. 0 divergences

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

write_csv(b_sum_conv, "outputs/trial_submodel/sum_rhat_ess.csv")

# Trace plots - global parameters (e.g. shrinkage, random effects)
trace_global <- traceplot(
  a_fit, 
  pars = c(
    "lp__", "mean_rate", "sd_rate",
    "sigma_a", "a_mean", "sigma_mu1", "mu1_mean",
    "tau", "shrinkage_effective_params"
  )
)

ggsave(
  "outputs/trial_submodel/trace_global.png", 
  trace_global, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots - latent mean for beta shrinkage
trace_z <- traceplot(a_fit, pars = "z")

ggsave(
  "outputs/trial_submodel/trace_shrinkage_coef.png", 
  trace_z, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots - local shrinkage
trace_lambda <- traceplot(a_fit, pars = "lambda")

ggsave(
  "outputs/trial_submodel/trace_shrinkage_lambda.png", 
  trace_lambda, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots - betas
trace_beta <- traceplot(a_fit, pars = c("beta_arm", "beta"))

ggsave(
  "outputs/trial_submodel/trace_beta.png", 
  trace_beta, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots - first 50 of primary parameters
trace_a <- traceplot(a_fit, pars = paste0("a[", 1:56, "]"))
trace_b <- traceplot(a_fit, pars = paste0("b[", 1:56, "]"))

ggsave(
  "outputs/trial_submodel/trace_a.png", 
  trace_a, 
  width = 15, 
  height = 10,
  units = "in"
)

ggsave(
  "outputs/trial_submodel/trace_b.png", 
  trace_b, 
  width = 15, 
  height = 10,
  units = "in"
)

# Trace plots - baseline log-hazard rate (mu1)
trace_mu1 <- traceplot(a_fit, pars = paste0("mu1[", 1:56, "]"))

ggsave(
  "outputs/trial_submodel/trace_mu1.png", 
  trace_mu1, 
  width = 15, 
  height = 10,
  units = "in"
)

# Posterior checks -------------------------------------------------------------

# Compare distribution of attrition counts
ppc_dens <- ppc_dens_overlay(a_imp_df$r, as.matrix(a_fit, pars = "r_new")[1:361, ])
ppc_int <- ppc_intervals(a_imp_df$r, as.matrix(a_fit, pars = "r_new")[1:361, ])

ggsave(
  "outputs/trial_submodel/ppc_density.jpg", 
  ppc_dens, 
  width = 8, 
  height = 8,
  units = "in"
)

ggsave(
  "outputs/trial_submodel/ppc_intervals.jpg", 
  ppc_int, 
  width = 8, 
  height = 8,
  units = "in"
)

# Residuals
ppc_ppe_mean <- ppc_stat(
  y = a_imp_df$r, 
  yrep = as.matrix(a_fit, pars = "r_new")[1:361, ], 
  stat = "mean"
)

ppc_ppe_sd <- ppc_stat(
  y = a_imp_df$r, 
  yrep = as.matrix(a_fit, pars = "r_new")[1:361, ], 
  stat = "sd"
)

ggsave(
  "outputs/trial_submodel/ppc_pred_err_mean.jpg", 
  ppc_ppe_mean, 
  width = 8, 
  height = 8,
  units = "in"
)

ggsave(
  "outputs/trial_submodel/ppc_pred_err_sd.jpg", 
  ppc_ppe_sd, 
  width = 8, 
  height = 8,
  units = "in"
)

# Pairwise correlation - coefficients
ppc_cor_beta <- a_draws[, grepl("^beta", names(a_draws))]
z_order_beta <- c("beta_arm", paste0("beta_", 1:30))

ppc_cor_beta_df <- as.data.frame(cor(ppc_cor_beta)) %>% 
  rownames_to_column(var = "row") %>% 
  pivot_longer(-row, names_to = "column", values_to = "cor") %>% 
  mutate(
    across(c(row, column), ~ gsub("\\[", "_", .)),
    across(c(row, column), ~ gsub("\\]", "", .)),
    across(c(row, column), ~ factor(., levels = z_order_beta))
  )

ppc_cor_beta_plot <- ppc_cor_df %>% 
  ggplot(aes(x = row, y = column, fill = cor)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Beta rows", y = "Beta columns", fill = "R")

ggsave(
  "outputs/ppc_cor_beta.jpg",
  ppc_cor_beta_plot,
  width = 8,
  height = 5,
  unit = "in"
)

# Shrinkage summary - 15 active covariates
hist(extract(a_fit, "shrinkage_effective_params")$shrinkage_effective_params)

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
  path <- paste0("outputs/trial_submodel/", save_name)
  plot <- ppc_plots[[name]]
  ggsave(path, plot, width = 6, height = 4, unit = "in")
  print(plot)
}
