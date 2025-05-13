
source("mars/design_model_mk2/scripts/00_config_functions_packages.R")

# Import data, model fit and posterior summaries
a_imp_df <- readRDS("mars/design_model_mk2/stan_dat.rds")
a_imp_mdl <- readRDS("mars/design_model_mk2/outputs/design2_fit.rds")
a_imp_sum <- readRDS("mars/design_model_mk2/outputs/design2_sum_draws.rds")

# Get beta names
b_beta_nm <- colnames(a_imp_df$X)
b_beta_nm <- c("trt_arm", b_beta_nm)

# Evaluate predictor significance for betas
b_betas <- a_imp_sum %>% 
  filter(grepl("^beta\\[|beta_arm", variable)) %>% 
  mutate(beta_name = factor(b_beta_nm, levels = b_beta_nm))

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
  "outputs/stan_validation.png",
  plot_b_betas,
  width = 8,
  height = 8,
  units = "in"
)
