
# Pool variables across multiarm trials
source("scripts/00_config.R")
library(tidyverse)

# Import bl chars
a_imp_bl <- readRDS("processed_data/agg_bl_fourth_ext.rds")
a_imp_set <- readRDS("processed_data/latest_dataset.rds")

# Add class
a_imp_class <- readRDS("processed_data/agg_bl_first_ext.rds")

a_add_class <- a_imp_class %>% 
  select(trial_id:trtcls5) %>% 
  left_join(a_imp_bl)

# Manually add baseline chars for two trials not in set
# write_csv(a_add_class, "scratch_data/add_n2_manually.csv")
a_imp_fix <- a_add_class %>% 
  filter(!trial_id %in% c("NCT02730377", "NCT02911948")) %>% 
  full_join(read_csv("scratch_data/added_n2_manually.csv")) %>% 
  arrange(trial_id, arm_id)

# Add treatment type
# Remove experimental arms
# Collapse trials with placebo and ac into one trttype
b_add_trttype <- a_imp_fix %>% 
  mutate(
    trttype = case_when(
      class %in% c("dpp4", "glp1", "sglt2") ~ NA,
      class == "placebo" ~ 0L,
      !class == "placebo" & !class %in% c("dpp4", "glp1", "sglt") ~ 1L
    )
  ) %>% 
  filter(!is.na(trttype)) %>% 
  group_by(trial_id) %>% 
  mutate(
    trttype = case_when(
      length(unique(trttype)) > 1 ~ 2L,
      TRUE ~ trttype
    )
  ) %>% 
  select(trial_id:arm_id, trttype, everything())

saveRDS(b_add_trttype, "scratch_data/ryan_bl_dataset.rds")

# Helper function - compute pooled SD when summarising group chars
pooledSD <- function(mu, delta, n) {
  # Get weighted mean
  keep_mu <- is.finite(mu) & is.finite(n)
  mu <- mu[keep_mu]
  n_mu <- n[keep_mu]
  if (length(n_mu) == 0) return(NA_real_)
  N <- sum(n_mu)
  if (N <= 1) return(NA_real_)
  Mu <- sum(n_mu * mu) / N
  
  # Get pooled SD
  keep_delta <- is.finite(delta) & is.finite(n)
  within_sd <- sum( (n[keep_delta] - 1) * (delta[keep_delta]^2) )
  bw_sd <- sum( n_mu * (mu - Mu)^2 )
  sqrt( (within_sd + bw_sd) / (N - 1) )
}

# Add sample size to get weighted means and sds collapsing trttype
b_imp_n <- b_add_trttype %>% 
  left_join(read_csv("scratch_data/aact_n.csv")) %>% 
  select(trial_id:trtcls5, n, everything())

# Manually add two missing
# write_csv(
#   b_imp_n %>% filter(is.na(n)),
#   "scratch_data/aact_n_manual.csv"
# )

# Add manual n
b_add_n2 <- b_imp_n %>% 
  filter(!trial_id %in% c("NCT02730377", "NCT02911948")) %>% 
  full_join(read_csv("scratch_data/aact_n_manual_added.csv")) %>% 
  arrange(trial_id, arm_id)

# Summarise across trttype per trial
# Weighted mean and sd accounting for within-group variation
# Sum race counts
b_sum_trttype <- b_add_n2 %>% 
  group_by(trial_id, trttype) %>% 
  reframe(
    arm_n = sum(n),
    across(asian:white, ~ weighted.mean(., n)),
    male = weighted.mean(male),
    across(current:former_never, ~ weighted.mean(., )),
    across(contains("mu"), ~ weighted.mean(., n)),
    age_d = pooledSD(age_mu, age_d, n),
    sbp_d = pooledSD(sbp_mu, sbp_d, n),
    bmi_d = pooledSD(bmi_mu, bmi_d, n),
    t2dur_d = pooledSD(t2dur_mu, t2dur_d, n),
    egfr_d = pooledSD(egfr_mu, egfr_d, n),
    hba1c_d = pooledSD(hba1c_mu, hba1c_d, n),
    hdl_d = pooledSD(hdl_mu, hdl_d, n),
    total_chol_d = pooledSD(total_chol_mu, total_chol_d, n)
  ) %>% 
  mutate(across(age_mu:total_chol_d, ~ round(., 2))) %>% 
  select(
    trial_id:trttype, arm_n, contains("age"), contains("sbp"), contains("bmi"),
    contains("t2dur"), contains("egfr"), contains("hba1c"), contains("hdl"),
    contains("total_chol"), everything()
  )

# Join remaining groups from design dataset
# Fix one hba1c where sd was interpreted as an entry error (does look like it,
# but it's a percentage seemingly, repeated in other table)
# Setting SD to NA, not clear what it could be
b_add_vars <- b_sum_trttype %>% 
  arrange(trial_id) %>% 
  cbind(
    a_imp_set %>% 
      arrange(trial_id) %>% 
      select(-trial_id, -n, -arm_n, -trttype, -attr_prop)
  ) %>% 
  as_tibble() %>% 
  select(source, trial_id, trttype, attr, everything()) %>% 
  mutate(hba1c_d = if_else(trial_id == "NCT01341717", NA_real_, hba1c_d))

# Check missingness
b_add_vars %>% 
  reframe(across(trial_id:former_never, ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  mutate(prop = value/344) %>% 
  arrange(desc(prop)) %>% 
  print(n = 25)

# Final checks for race - manual look, likely some of these trials have
# homogeneous populations
# b_add_vars %>% 
#   filter(if_all(c(asian, black, other, white), ~ is.na(.))) %>% 
#   write_csv("scratch_data/race_check_homogen.csv")

# Add in fixes - most were homogeneous populations, mostly asian
b_imp_fix <- read_csv("scratch_data/race_check_homogen.csv")

b_fix_race <- b_add_vars %>% 
  filter(!trial_id %in% b_imp_fix$trial_id) %>% 
  full_join(b_imp_fix) %>% 
  mutate(across(age_mu:former_never, ~ round(., 2))) %>% 
  select(-c(
    age_lower_dbl:age_upper_dbl, 
    bmi_lower:bmi_upper,
    male_prop:egfr_mean
  )) %>% 
  mutate(
    white = case_when(
      asian + black + other + white == 1.01 ~ white - 0.01, 
      asian + black + other + white == 0.99 ~ white + 0.01, 
      TRUE ~ white
    )
  ) %>% 
  select(source:former_never, everything())

# Recheck missingness
b_fix_race %>% 
  reframe(across(trial_id:former_never, ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  mutate(prop = value/344) %>% 
  arrange(desc(prop)) %>% 
  print(n = 25)

# Save
saveRDS(b_fix_race, "processed_data/agg_bl_pooled_arms.rds")
