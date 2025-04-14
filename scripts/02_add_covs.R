
# Extract design features as covariates from aact datasets
source("scripts/00_functions_packages_etc.R")

a_imp_df <- read_csv("processed_data/df_counts_arms.csv")
a_imp_covs <- read_csv("data/aact_covars.csv")

# Get design features for trials in set
b_join_covs <- a_imp_df %>% left_join(a_imp_covs %>% rename(trial_id = nct_id))
b_join_covs %>% count(trial_id, source) %>% count(source)

# Summarise missingness (change missing count when all included)
# b_join_covs %>% 
#   reframe(across(everything(), ~ sum(is.na(.)))) %>% 
#   pivot_longer(
#     everything(),
#     names_to = "var",
#     values_to = "na"
#   ) %>% 
#   mutate(missing = if_else(na == 346, 1, 0)) %>% 
#   write_csv("scratch_data/variable_missingness.csv")

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
  ) %>% 
  rename(trttype = type)

b_imp_manual %>% count(trial_id, source) %>% count(source)

# Add new trials
# Add n for new ipd trials
z_n <- b_add_new %>% filter(is.na(n)) %>% select(trial_id)
z_n$n_new <- c(361, 680, 484, 566, 243, 448, 275)

b_add_new <- b_imp_manual %>% 
  full_join(b_join_covs %>% filter(!trial_id %in% b_imp_manual$trial_id)) %>% 
  filter(trial_id %in% b_join_covs$trial_id) %>% 
  left_join(z_n, by = "trial_id") %>% 
  mutate(n = coalesce(n, n_new)) %>% 
  select(-n_new)

b_add_new %>% count(trial_id, source) %>% count(source)
b_add_new %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save version with new trials added to manually find missing data
# write_csv(b_add_new, "scratch_data/aact_manual_find_missingness_new.csv")

# No changes needed, only the 16 with no facilities are missing

# Impute missing facilities based on the sample size
library(mice)

c_impute_vars <- b_add_new %>% 
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
  left_join(b_add_new %>% select(-n, -n_facilities)) %>% 
  select(trial_id, source:number_of_arms, n:n_facilities, everything())

# Final missingness check
c_impute_join %>% reframe(across(everything(), ~ sum(is.na(.))))
c_impute_join %>% count(trial_id, source) %>% count(source)

# Remove unused variables
c_impute_join <- c_impute_join %>% select(-c(n_soms:outcomes_assessor_masked))

# Save
# write_csv(c_impute_join, "processed_data/df_counts_arms_covs.csv")

# Format for modelling
b_fmt <- c_impute_join %>% 
  mutate(
    phase = if_else(!grepl("Phase", phase), paste0("Phase ", phase), phase)
  ) %>% 
  pivot_wider(names_from = "phase", values_from = "phase") %>% 
  pivot_wider(names_from = "masking", values_from = "masking") %>% 
  mutate(
    across(`Phase 3`:`Phase Not Applicable`, ~ if_else(is.na(.), 0, 1)),
    trttype = if_else(trttype == "placebo", 1, 0),
    across(Double:Triple, ~ if_else(is.na(.), 0, 1))
  ) %>% 
  rename_at(vars(`Phase 3`:`Phase Not Applicable`), ~ gsub(" ", "_", .)) %>% 
  rename_at(vars(Double:Triple), ~ gsub(" ", "_", .)) %>% 
  rename_at(vars(Phase_3:Triple), ~ tolower(.)) %>% 
  rename(
    n_arms = number_of_arms,
    phase_notapp = phase_not_applicable,
    open_label = `none_(open_label)`
  ) %>% 
  select(
    trial_id:source, n, arm_n:n_arms, phase_3:phase_notapp,
    double:triple, everything()
  ) %>% 
  mutate(
    across(
      c(phase_3:triple, region_of_the_americas:eastern_mediterranean_region),
      ~ as.integer(.)
    ),
    across(
      region_of_the_americas:eastern_mediterranean_region,
      ~ if_else(. > 1, 1, .)
    )
  )

# Final checks
b_fmt %>% count(trial_id, source) %>% count(source)
b_fmt %>% reframe(across(everything(), ~ sum(is.na(.)))) %>% View()

# Save
write_csv(b_fmt, "processed_data/df_stan_n348.csv")
