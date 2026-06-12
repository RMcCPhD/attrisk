
source("scripts/00_config.R")
library(tidyverse)

# Import trial covariates
a_imp_df <- read_csv("processed_data/definitive_dataset.csv")

# Get groups (347 trials)
a_grps <- a_imp_df %>% select(trial_id) %>% distinct()

# Import eligibility criteria
a_imp_elig <- read_csv("data_old/elig.csv")

# Get elig for 347 trials
b_elig_trials <- a_imp_elig %>% inner_join(a_grps)
b_elig_trials %>% count(trial_id)
glimpse(b_elig_trials)

# Find one trial not present
# NCT04170998 - add manually at end, some on ctgov
a_grps %>% anti_join(b_elig_trials %>% select(trial_id))

# Take age, heart failure, comorbidities, con meds and disability flags
b_elig_ext <- b_elig_trials %>% 
  select(
    trial_id, minimum_age:maximum_age, 
    alcohol_and_drugs:other_medicines
  )

# Missingness and variable selection
b_elig_ext %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "nas") %>% 
  mutate(nas_pcnt = nas/350)

# Harmonise variables
c_harmonise <- b_elig_ext %>% 
  mutate(
    # Heart failure (binary - include/exclude) - DNE as reference
    hf = case_when(
      grepl("Do not exclude", heart_failure) ~ 0L,
      grepl("Class", heart_failure) ~ 1L
    )
  ) %>% 
  select(trial_id, contains("age"), hf)

# Import eligibility limits for bmi, bp, durn and egfr
d_imp_files <- list.files("data", pattern = "upper.csv$")
d_imp_paths <- list.files("data", pattern = "upper.csv$", full.names = TRUE)
d_imp_lst <- lapply(d_imp_paths, read.csv)
names(d_imp_lst) <- d_imp_files

# Tidy bmi
d_bmi <- d_imp_lst$bmi_lower_upper.csv %>% 
  rename(trial_id = nct_id) %>% 
  inner_join(a_grps) %>% 
  select(trial_id, inex, operator, lower, upper) %>% 
  separate(operator, into = c("op_lower", "op_upper"), sep = ";") %>% 
  mutate(
    op_upper = if_else(is.na(op_upper) & !is.na(op_lower), op_lower, op_upper),
    lower = paste0(op_lower, lower),
    upper = paste0(op_upper, upper),
    across(lower:upper, ~ if_else(grepl("NA", .), NA, .))
  ) %>% 
  select(-contains("op_"))

# Get inclusion for exclusion statements in bmi
d_bmi_all_incl <- d_bmi %>% 
  as_tibble() %>% 
  mutate(
    lower_new = case_when(
      inex == "exclusion" & grepl("<|<=", lower) ~ gsub("<|<=", ">", lower),
      TRUE ~ lower
    ),
    upper_new = case_when(
      inex == "exclusion" & grepl(">|>=", upper) ~ gsub(">|>=", "<", upper),
      TRUE ~ upper
    )
  ) %>% 
  select(-c(inex, lower, upper)) %>% 
  rename(
    bmi_lower = lower_new,
    bmi_upper = upper_new
  )

# Join bmi
d_bmi_join <- c_harmonise %>% left_join(d_bmi_all_incl)

# Tidy egfr
e_egfr <- d_imp_lst$egfr_lower_upper.csv %>% 
  rename(trial_id = nct_id) %>% 
  inner_join(a_grps) %>% 
  select(trial_id, lower, upper) %>% 
  rename(
    egfr_lower = lower,
    egfr_upper = upper
  )

# Join egfr
e_egfr_join <- d_bmi_join %>% 
  left_join(e_egfr) %>% 
  rename(
    age_lower = minimum_age,
    age_upper = maximum_age
  ) %>% 
  select(trial_id, hf, everything())

# Save
write_csv(e_egfr_join, "processed_data/elig_tidy.csv")
