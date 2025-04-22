
source("scripts/00_functions_packages_etc.R")

# Import posterior summaries and stan dataset
a_imp_trial_sum <- readRDS("outputs/posterior_sum_trials.rds")
a_imp_global_sum <- readRDS("outputs/posterior_sum_global.rds")
a_imp_df <- readRDS("processed_data/stan_dat.rds")

# Get covariate names
z_cov_names <- colnames(a_imp_df$X)

# Join names to identify coefficients
b_betas <- a_imp_global_sum %>% 
  filter(grepl("beta", variable)) %>% 
  select(variable, mean) %>% 
  mutate(n = if_else(variable == "beta_arm", 0, row_number() - 1)) %>% 
  left_join(
    z_cov_names %>% 
      as.data.frame() %>% 
      rename(cov_name = ".") %>% 
      mutate(n = row_number())
  ) %>% 
  mutate(cov_name = if_else(is.na(cov_name), "trttype", cov_name))

saveRDS(b_betas, "outputs/betas_covs.rds")
