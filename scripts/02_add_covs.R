
# Extract design features as covariates from aact datasets
source("scripts/00_functions_packages_etc.R")

a_imp_df <- read_csv("processed_data/df_counts_arms.csv")
a_imp_covs <- read_csv("data/aact_covars.csv")

# TEMP remove IPD NAs (n=346)
a_imp_df <- a_imp_df %>% filter(!is.na(arm_n_attr))

# Get design features for trials in set
b_join_covs <- a_imp_df %>% left_join(a_imp_covs %>% rename(trial_id = nct_id))

# Summarise missingness (change missing count when all included)
b_join_covs %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(
    everything(),
    names_to = "var",
    values_to = "na"
  ) %>% 
  mutate(missing = if_else(na == 346, 1, 0)) %>% 
  write_csv("scratch_data/variable_missingness.csv")

# Save joined design features to fill missingness manually
# write_csv(b_join_covs, "scratch_data/aact_manual_find_missingness.csv")

# Import with manually imputed design features
# Found on ctgov and published materials
# Set individual regions as 1 or 0 to simplify variable
# Number of regions is primary variable for this context
b_imp_manual <- read_csv("scratch_data/aact_manual_found_missingness.csv") %>% 
  mutate(
    across(
      region_of_the_americas:eastern_mediterranean_region,
      ~ if_else(. >= 1, 1, 0)
    )
  )

# Impute missing facilities based on the sample size
library(mice)

c_impute_vars <- b_imp_manual %>% 
  select(trial_id, n, n_facilities) %>% 
  distinct() 

c_impute <- mice(
  c_impute_vars[, c("n", "n_facilities")],
  m = 5, 
  method = "pmm", 
  seed = 123
)

summary(c_impute)

c_impute_complete <- complete(c_impute, 1)
c_impute_complete <- c_impute_complete %>% cbind(trial_id = c_impute_vars$trial_id)

# Join imputed values
c_impute_join <- c_impute_complete %>% 
  left_join(b_imp_manual %>% select(-n, -n_facilities)) %>% 
  select(trial_id, source:number_of_arms, n:n_facilities, everything())

# Final missingness check
c_impute_join %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save
write_csv(c_impute_join, "processed_data/df_counts_arms_covs.csv")
