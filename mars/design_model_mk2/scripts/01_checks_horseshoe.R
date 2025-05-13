
source("mars/design_model_mk2/scripts/00_config_functions_packages.R")

# Import data, model fit and posterior summaries
a_imp_df <- readRDS("mars/design_model_mk2/stan_dat.rds")
a_imp_mdl <- readRDS("mars/design_model_mk2/outputs/horseshoe/design2_fit.rds")
a_imp_sum <- readRDS("mars/design_model_mk2/outputs/horseshoe/design2_sum_draws.rds")

# Compute WAIC for model
a_draws <- as.matrix(a_imp_mdl, pars = "log_lik")
loo(a_draws)
loo::waic(a_draws)

# Get beta names
b_beta_nm <- colnames(a_imp_df$X)
b_beta_nm <- c("trt_arm", b_beta_nm)

# Evaluate predictor significance for betas
b_betas <- a_imp_sum %>% 
  filter(grepl("^beta\\[|beta_arm", variable)) %>% 
  mutate(beta_name = factor(b_beta_nm, levels = b_beta_nm))

# Import and join betas with shrinkage
c_both <- b_betas %>% 
  mutate(shrinkage = "horseshoe") %>% 
  full_join(
    read_csv("mars/design_model_mk2/scratch_shrinkage.csv") %>% 
      mutate(shrinkage = "gaussian g-l")
  )

plot_b_betas <- c_both %>% 
  ggplot(
    aes(
      x = mean, xmin = `2.5%`, xmax = `97.5%`, 
      y = fct_rev(beta_name), 
      colour = shrinkage,
    )
  ) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0.05, colour = "orange", linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = -0.05, colour = "orange", linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(n.breaks = 20) +
  labs(x = "mean (95% credible intervals) effect", y = "covariate")

ggsave(
  "outputs/stan_comparison.png",
  plot_b_betas,
  width = 8,
  height = 8,
  units = "in"
)
