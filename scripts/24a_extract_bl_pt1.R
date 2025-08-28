
# First part of extracting aggregated baseline characteristics
# Takes in in tidied aggregated baseline characteristics from two repos
# of Type2DiabetesSystematicReview (additional_covariates, nma_more_covariates)

# Variables of interest
# Trttype, smoke, sex, sbp
# Race (asian, black/aa, white, other)
# Number of con meds and comorbidities
# HbA1c, eGFR, duration of T2DM
# Total and HDL cholesterol
# BMI, age10

# Obtained
# Race, age, sex, bmi, t2dur, smoke, egfr, sbp, chol

# To get
# Trttype, hba1c, hdl

source("scripts/00_config.R")
library(tidyverse)

# Function to check missingness
checkNA <- function(df) {
  df %>% 
    reframe(across(everything(), ~ sum(is.na(.)))) %>% 
    pivot_longer(everything()) %>% 
    filter(value > 0)
}

# Import trial set
a_imp_set <- readRDS("processed_data/latest_dataset.rds")

# Add arm ids
a_set_grps <- a_imp_set %>% 
  distinct(trial_id) %>% 
  left_join(read_csv("created_metadata/arm_id_trt.csv"))

# additional_covariates/Scratch_data - Race ------------------------------------

path_race <- "../additional_covariates/Scratch_data/aact_baseline_race.Rds"
b_imp_race <- readRDS(path_race)
checkNA(b_imp_race) # 0 NAs
b_imp_race_lkp <- read_csv("created_metadata/race_new.csv")

# Create lookup for simplifying race
# b_imp_race %>% count(var_level) %>% write_csv("scratch_data/race_create_lkp.csv")

b_race_tidy <- b_imp_race %>% 
  select(nct_id, arm_id, var_level, p) %>% 
  left_join(b_imp_race_lkp %>% select(-n)) %>% 
  select(nct_id:arm_id, race, p) %>% 
  group_by(nct_id, arm_id, race) %>% 
  reframe(p = sum(p)) %>% 
  pivot_wider(names_from = race, values_from = p)

# 345 trials
b_race_tidy %>% count(nct_id)

# !! 253/344 in attrisk set !!
b_race_tidy %>% 
  inner_join(
    a_imp_set %>% 
      distinct(trial_id) %>% 
      rename(nct_id = trial_id)
  ) %>% 
  count(nct_id)

# Join with trial set IDs
b_race_join <- a_set_grps %>% 
  left_join(b_race_tidy %>% rename(trial_id = nct_id))

# additional_covariates/Scratch_data - age, sex, bmi, t2dur --------------------

path1_age_sex_bmi_t2dur <- "../additional_covariates/Scratch_data/"
path2_age_sex_bmi_t2dur <- "age_sex_bmi_duration_agg_forhba1c_reg.rds"
path_age_sex_bmi_t2dur <- paste0(path1_age_sex_bmi_t2dur, path2_age_sex_bmi_t2dur)
b_age_sex_bmi_t2dur <- readRDS(path_age_sex_bmi_t2dur)
checkNA(b_age_sex_bmi_t2dur)

b_age_sex_bmi_t2dur_tidy <- bind_rows(b_age_sex_bmi_t2dur$agg_durn_bmi) %>% 
  mutate(male_p = male/n) %>% 
  select(
    nct_id, arm_id_unq, arm_lvl, trtcls5, age_m, age_sd, male = male_p,
    bmi_mean, bmi_sd, dur_mean, dur_sd
  )

# !! 224/344 in attrisk set !!
b_age_sex_bmi_t2dur_tidy %>% 
  inner_join(
    a_imp_set %>% 
      distinct(trial_id) %>% 
      rename(nct_id = trial_id)
  ) %>% 
  count(nct_id)

checkNA(b_age_sex_bmi_t2dur_tidy) # 0 NAs

b_age_sex_bmi_t2dur_join <- a_set_grps %>% 
  left_join(
    b_age_sex_bmi_t2dur_tidy %>% 
      rename(
        trial_id = nct_id,
        arm_id = arm_id_unq
      )
  )

# additional_covariates/Scratch_data - smoke -----------------------------------

b_smoke <- readRDS("../additional_covariates/Scratch_data/aact_smoke.rds")
checkNA(b_smoke)

# Error in one trial, all labelled current
# Fixing - current and former
b_smoke_fix <- b_smoke %>% 
  filter(nct_id == "NCT01528254") %>% 
  arrange(nct_id, arm_id) %>% 
  mutate(
    category_fix = c(0, 1, 0, 1, 0, 1),
    category_level = if_else(category_fix == 0, category_level, "former")
  )

# Seems to be an error for one trial - all labelled current
b_smoke_tidy <- b_smoke %>% 
  filter(!nct_id %in% b_smoke_fix$nct_id) %>% 
  full_join(b_smoke_fix) %>% 
  arrange(nct_id, arm_id) %>% 
  select(nct_id:p) %>% 
  mutate(
    category = case_when(
      grepl("former|never|unk", category_level) ~ "former_never",
      TRUE ~ "current"
    )
  ) %>% 
  group_by(nct_id, arm_id, category) %>% 
  reframe(p = sum(p)) %>% 
  pivot_wider(names_from = category, values_from = p)

# !! 12/344 in attrisk set !!
b_smoke_tidy %>% 
  inner_join(
    a_imp_set %>% 
      distinct(trial_id) %>% 
      rename(nct_id = trial_id)
  ) %>% 
  count(nct_id)

checkNA(b_smoke_tidy) # 19 NAs in current, 125 in former/never

b_smoke_join <- a_set_grps %>% 
  left_join(b_smoke_tidy %>% rename(trial_id = nct_id))

# additional_covariates/Scratch_data - egfr ------------------------------------

b_egfr <- readRDS("../additional_covariates/Scratch_data/egfr_params_m_s.rds")
checkNA(b_egfr) # No vital missingness

b_egfr_tidy <- b_egfr %>% 
  select(nct_id, arm_comb, variable, first, second) %>% 
  rename(mean = first, sd = second) %>% 
  group_by(nct_id, arm_comb) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = variable, values_from = c(mean, sd), names_vary = "slowest"
  )

# !! 50/344 in attrisk set !!
b_egfr_tidy %>% 
  inner_join(
    a_imp_set %>% 
      distinct(trial_id) %>% 
      rename(nct_id = trial_id)
  ) %>% 
  count(nct_id)

checkNA(b_egfr_tidy) # 0 NAs

b_egfr_join <- a_set_grps %>% 
  left_join(
    b_egfr_tidy %>% 
      rename(trial_id = nct_id) %>% 
      rename(arm_id = arm_comb)
  )

# additional_covariates/Scratch_data - sbp, egfr, chol -------------------------

path1_sbp_egfr_sbp_chol <- "../additional_covariates/Scratch_data/"
path2_sbp_egfr_sbp_chol <- "continuous_covariates.rds"
path_sbp_egfr_sbp_chol <- paste0(path1_sbp_egfr_sbp_chol, path2_sbp_egfr_sbp_chol)
b_sbp_egfr_sbp_chol <- readRDS(path_sbp_egfr_sbp_chol)
checkNA(b_sbp_egfr_sbp_chol) # 12 NAs in "second"

# Taking mean/SD for now
# Take first where there are 2 rows for a trial baseline char
b_sbp_chol_tidy <- b_sbp_egfr_sbp_chol %>% 
  filter(
    first_second == "mean_sd", 
    variable %in% c("blood_pressure", "hdl", "total_cholesterol")
  ) %>% 
  select(nct_id, arm_comb, variable, first, second) %>% 
  group_by(nct_id, arm_comb) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  rename(mean = first, sd = second) %>% 
  pivot_wider(
    names_from = variable, values_from = c(mean, sd), names_vary = "slowest"
  )

# !! 103/344 in attrisk set, not complete cases !!
b_sbp_chol_tidy %>% 
  inner_join(
    a_imp_set %>% 
      distinct(trial_id) %>% 
      rename(nct_id = trial_id)
  ) %>% 
  count(nct_id)

checkNA(b_sbp_chol_tidy) # 0 NAs

# None have hdl, removing
b_sbp_chol_join <- a_set_grps %>% 
  left_join(
    b_sbp_chol_tidy %>% 
      rename(trial_id = nct_id) %>% 
      rename(arm_id = arm_comb)
  ) %>% 
  select(-contains("hdl"))

# Combine extracts -------------------------------------------------------------

# Get joined datasets
lst_join <- lapply(ls(pattern = "_join"), function(x) get(x))

# Join together, tidy
join_df <- reduce(lst_join, left_join) %>% 
  group_by(class) %>% 
  fill(trtcls5, .direction = "updown") %>% 
  ungroup() %>% 
  select(-arm_lvl)

# Save
saveRDS(join_df, "processed_data/agg_bl_first_ext.rds")

