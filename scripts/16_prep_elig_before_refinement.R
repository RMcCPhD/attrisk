

source("scripts/00_config.R")
library(tidyverse)

# Import trial covariates
a_imp_df <- read_csv("processed_data/definitive_dataset.csv")

# Get groups (347 trials)
a_grps <- a_imp_df %>% select(trial_id) %>% distinct()

# Import eligibility criteria
a_imp_elig <- read_csv("data/elig.csv")
str(a_imp_elig)

# Get elig for 347 trials
b_elig_trials <- a_imp_elig %>% inner_join(a_grps)
b_elig_trials %>% count(trial_id)

# Find one trial not present
# NCT04170998 - add manually at end, some on ctgov
a_grps %>% anti_join(b_elig_trials %>% select(trial_id))

# Missingness and variable selection
b_sum_nas <- b_elig_trials %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "nas") %>% 
  mutate(nas_pcnt = nas/350)

# write_csv(b_sum_nas, "scratch_data/elig_select_vars.csv")

# Get selected variables (based on missingness)
c_selected <- read_csv("scratch_data/elig_select_vars.csv") %>% 
  filter(flag == 1) %>% 
  select(name)

c_selected <- c_selected$name

# Convert t2dm all to months
c_select_vars <- b_elig_trials %>% 
  select(all_of(c_selected)) %>% 
  mutate(
    t2d_conv = as.numeric(gsub("<|>", "", t2dm_duration)),
    t2d_conv = case_when(
      t2dm_duration_unit == "days" ~ t2d_conv / 30,
      t2dm_duration_unit == "weeks" ~ t2d_conv / 4,
      t2dm_duration_unit == "years" ~ t2d_conv * 12,
      TRUE ~ t2d_conv
    ),
    t2d_conv = round(t2d_conv)
  ) %>% 
  select(trial_id:hba1c_max, contains("t2d"), everything())

# Function to quickly check vars for harmonising
check_var <- function(var) {
  freq <- c_select_vars %>% count({{ var }})
  print(freq)
}

# Harmonise variables
d_harmonise <- c_select_vars %>% 
  mutate(
    # Min age (binary) - 18 as reference
    min_age = case_when(
      minimum_age == 18 ~ 0L,
      minimum_age != 18 & !is.na(minimum_age) ~ 1L
    ),
    # Max age (binary) - >75 as reference
    max_age = case_when(
      maximum_age > 75 ~ 1L,
      maximum_age <= 75 ~ 0L
    ),
    # Min bmi (binary) - >75 as reference
    min_bmi = case_when(
      bmi > 75 ~ 1L,
      bmi <= 75 ~ 0L
    ),
    # Max bmi (binary) - >75 as reference
    max_age = case_when(
      bmi > 75 ~ 1L,
      bmi <= 75 ~ 0L
    ),
    # Min hba1c (binary) - >=7 as reference
    min_hba1c = case_when(
      hba1c_min >= 7 ~ 0L,
      hba1c_min < 7 ~ 1L
    ),
    # Max hba1c (binary) - >=10 as reference
    max_hba1c = case_when(
      hba1c_max >=10 ~ 0L,
      hba1c_max < 10 ~ 1L
    ),
    # T2D duration (binary) - <6m as reference
    t2d_lim = case_when(
      t2d_conv < 6 ~ 0L,
      t2d_conv >= 6 ~ 1L
    ),
    # Alcohol and drugs (binary) - No as reference
    alc_drg = case_when(
      alcohol_and_drugs %in% c("No", "no") ~ 0L,
      alcohol_and_drugs %in% c("Yes", "yes") ~ 1L
    ),
    # Heart failure (binary) - DNE as reference
    hf = case_when(
      grepl("Do not exclude", heart_failure) ~ 0L,
      grepl("Class", heart_failure) ~ 1L
    ),
    # Comorbidities (binary) - yes as reference
    comorb = case_when(
      !grepl("none", non_cardiovascular_comorbidities) ~ 0L,
      TRUE ~ 1L
    ),
    # Medications (binary) - yes as reference
    meds = case_when(
      !grepl("none", other_medicines) ~ 0L,
      TRUE ~ 1L
    ),
    # Physical disability - no as reference
    phys = case_when(
      grepl("No|no", physical_disability) ~ 0L,
      TRUE ~ 1L
    ),
    # Mental disability - no as reference
    mh = case_when(
      grepl("No|no", mental_disability) ~ 0L,
      TRUE ~ 1L
    )
  )

# Select variables
e_vars <- d_harmonise %>% 
  select(trial_id, min_age:mh) %>% 
  mutate(across(-trial_id, as.integer))

# Save (manually address NAs where possible and add missing trial)
write_csv(e_vars, "processed_data/elig_tidy.csv")

# Join with rest of covariates after adding NAs/missing trial