
source("scripts/00_functions_packages_etc.R")

# 348 trials, 86 have fu ratios (IPD trials)
a_imp_df <- read_csv("processed_data/df_stan_n348.csv")
a_imp_params <- read_csv("data/parameters.csv")

# Get IPD trials
b_ipd <- a_imp_df %>% filter(source == "ipd")

# Get Gompertz parameters (33 trials)
b_params <- a_imp_params %>% 
  filter(dist == "gompertz") %>% 
  select(-se) %>% 
  pivot_wider(names_from = "parameter", values_from = "est") %>% 
  rename(trial_id = ctgov) %>% 
  filter(trial_id %in% b_ipd$trial_id)



