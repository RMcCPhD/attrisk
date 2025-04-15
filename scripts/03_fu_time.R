
# Add mean follow-up time
# Available for 86 IPD trials, approximated for remaining
source("scripts/00_functions_packages_etc.R")

# 348 trials, 86 have fu ratios (IPD trials)
a_imp_df <- read_csv("processed_data/df_stan_n348.csv")
a_imp_fu <- read_csv("data/fu_ratio_n86.csv")

# Prepare fu ratio
b_fu <- a_imp_fu %>% 
  separate(
    collapsed,
    into = c("trttype", "mft_ratio"),
    sep = " "
  )
