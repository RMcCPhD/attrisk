
# Prepare eligibility criteria that were incorporated into development
# i.e. Included as trial design features

source("scripts/00_config.R")
library(tidyverse)

# Import harmonised eligibility criteria
a_imp_elig <- read_csv("data/elig.csv")

# Get AgD trial set
a_imp_agd <- readRDS("processed_data/latest_dataset.rds")

# Check how many AgD trials are in eligibility (should be all)
# 1 is not present (NCT04170998)
a_imp_agd %>% anti_join(a_imp_elig %>% distinct(trial_id))

# Join and isolate trial that is not included
# Select relevant variables
b_elig_agd <- a_imp_agd %>% 
  distinct(trial_id) %>% 
  inner_join(a_imp_elig) %>% 
  select(
    trial_id, 
    contains(c("_age", "bmi", "renal"))
  )

# Manually create row for NCT04170998
# Had added these in for model development
b_missing <- b_elig_agd[1,]
# write_csv(b_missing, "scratch_data/elig_missing_trial.csv")

# Full-join NCT04170998
b_elig_all <- b_elig_agd %>%
  mutate(across(everything(), as.character)) %>% 
  full_join(
    read_csv("scratch_data/elig_missing_trial.csv") %>% 
      mutate(across(everything(), as.character))
  ) %>% 
  arrange(trial_id) %>% 
  rename(nct_id = trial_id)

# Save
write_csv(b_elig_all, "../../Attrisk/GitHub/eligibility_criteria.csv")