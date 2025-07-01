
# Imputing missing participant characteristics
# Previously imputed age, bmi, t2d duration, hba1c and male proportion
# This was insufficient - only used arm_n as a covariate, and invalid to impute
# male proportion
# Start fresh in this script for the participant characteristics

source("scripts/00_config.R")
library(tidyverse)

# Import latest version of aggregated dataset
a_imp_df <- read_csv("processed_data/definitive_dataset_with_elig.csv")

# Remove participant characteristics
# 347 trials
b_rm_pct <- a_imp_df %>% select(-c(scale_male_prop:hba1c))

# Import latest version of covariates from t2dSR repo
# Get attrisk trials
# 345 trials have data
c_imp_covs <- read_csv("data/base_dsp.csv") %>% 
  inner_join(
    b_rm_pct %>% 
      select(trial_id) %>% 
      distinct()
  )

# Missing trial overall (differs by covariate) - NCT02730377
b_rm_pct %>% anti_join(c_imp_covs %>% select(trial_id) %>% distinct())

# Tidy characteristics ---------------------------------------------------------

# Take relevant variables
d_vars <- c_imp_covs %>% 
  select(trial_id:arm_id, variable, measure_type:second, category_level)

# 2 age missing measure_type, but age always mean where non-NA
d_vars %>% count(variable)
d_vars %>% count(variable, measure_type)

# 2 age where second measure is labelled "median"
# Majority of male prop is a count, some only percentage
d_vars %>% count(variable, measure_type, first_format, second_format)

# Dealing with each covariate individually, then joining together
# Ignoring ethnicity, two trials have it

# Sample size  - 1 trial missing
e_sample <- d_vars %>% 
  filter(variable == "n") %>% 
  select(trial_id:arm_id, sample = first)

# Male - 7 trials missing
# Two missing percentage, keep n to identify
# Calculate n where missing by joining sample sizes
# Keep percentage in case missing n can be identified where sample is missing
e_male <- d_vars %>% 
  filter(variable == "male") %>% 
  select(trial_id:arm_id, first, second) %>% 
  rename(
    male_n = first,
    male_pcnt = second
  ) %>% 
  mutate(male_n = if_else(!is.na(male_n) & is.na(male_pcnt), male_n, NA)) %>% 
  left_join(e_sample) %>% 
  mutate(
    male_n = case_when(
      is.na(male_n) & !is.na(sample) ~ round(sample * (male_pcnt/100)),
      TRUE ~ male_n
    )
  ) %>% 
  select(-sample)

# Age - 9 trials missing
# Set median for second measure as NA (not usable, rest are sd)
e_age <- d_vars %>% 
  filter(variable == "age") %>% 
  mutate(
    second = case_when(
      second_format == "median" ~ NA,
      TRUE ~ second
    )
  ) %>% 
  rename(
    age_mean = first,
    age_sd = second
  ) %>% 
  select(trial_id:arm_id, age_mean, age_sd)

# Race - 119 trials missing
# Import lookup from t2d sr repo
e_race_lkp <- read_csv("created_metadata/race_smplfy.csv")

# Calculate n using sample and race percentages
# Simplify names and remove collapsed races (name already harmonised by lookup)
e_race <- d_vars %>% 
  filter(variable == "race") %>% 
  select(trial_id:arm_id, category_level, second) %>% 
  rename(
    race = category_level,
    race_pcnt = second
  ) %>% 
  mutate(race_pcnt = if_else(race_pcnt > 100, 100, round(race_pcnt, 2))) %>% 
  left_join(e_sample) %>% 
  mutate(n = case_when(!is.na(sample) ~ round(sample * (race_pcnt/100)))) %>% 
  left_join(e_race_lkp %>% rename(race = category_level)) %>% 
  select(-sample) %>% 
  filter(!grepl("collapsed", race)) %>% 
  rename(race_hrm = smplfy)

# Join together
f_join <- e_sample %>% 
  left_join(e_male) %>% 
  left_join(e_age) %>% 
  left_join(
    e_race %>% 
      select(-c(race, race_pcnt)) %>% 
      mutate(
        race_hrm = case_when(
          race_hrm %in% c("a_ai_h_pi", "hispanic_latino", "other", "unknown") ~ "other",
          TRUE ~ race_hrm
        )
      ) %>% 
      group_by(trial_id, arm_id, race_hrm) %>% 
      reframe(n = sum(n)) %>% 
      pivot_wider(names_from = "race_hrm", values_from = "n")
  )

# Most missingness in race
f_join %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save before imputation
write_csv(f_join, "processed_data/participant_chars_before_imputation.csv")

# Impute missing continuous variables ------------------------------------------