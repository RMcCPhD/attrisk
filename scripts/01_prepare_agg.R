
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## SUMMARY OF SCRIPT

#' ## 1. PREPARE AACT DATA - OVERALL ATTRITION
#' Get aact trials
#' Get overall attrition counts (ignoring arms for first initial stan model)
#' Add design variables
#' Check missingness, remove variables with high missingness

#' ## 2. PREPARE AACT DATA - ARM-LEVEL ATTRITION
#' TBC

#' END RESULT
#' 207 aact trials with 12 design variables (excl. duration)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

source("scripts/00_functions_packages_etc.R")

## 1. GET ELIGIBLE TRIALS FOR PROJECT ##########################################

# Import data - attrition counts, design variables
imp_trials <- read_rds("data/trials_n392_attrisk.rds")
imp_design <- read_rds("data/design_vars_n372.rds")

# 300 agg, 92 ipd
imp_trials %>% count(trial_id, source) %>%  count(source)

# 7 ipd trials missing arm attrition counts, removing for test
imp_trials %>% group_by(trial_id) %>% filter(any(is.na(arm_n_attr)))

# 385 trials
imp_trials_counts <- imp_trials %>% 
  group_by(trial_id) %>% 
  filter(!any(is.na(arm_n_attr))) %>% 
  ungroup() %>% 
  select(-c(has_ipd:on_aact))

# Add overall attrition counts
add_oattr <- imp_trials_counts %>% 
  group_by(trial_id) %>% 
  mutate(attr = sum(arm_n_attr)) %>% 
  ungroup()

# Get overall attrition dataset with design variables
df_oattr <- add_oattr %>% 
  select(trial_id, n, attr) %>% 
  distinct() %>% 
  rename(nct_id = trial_id) %>% 
  left_join(imp_design)

# Save without any removals for missingness
saveRDS(df_oattr, "processed_data/df_oattr_all_design_vars.rds")

# Check missingness
remove_vars <- df_oattr %>% 
  reframe(across(-c(nct_id:attr), ~ sum(is.na(.)))) %>% 
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "na"
  ) %>% 
  arrange(desc(na)) %>% 
  filter(na > 2)

remove_vars <- remove_vars$variable

# Remove variables with high missingness
df_oattr_rem_vars <- df_oattr %>% select(-all_of(remove_vars))

# Save
saveRDS(df_oattr_rem_vars, "processed_data/df_oattr.rds")




