
# Prepare eligibility criteria that were incorporated into development
# i.e. Included as trial design features

source("scripts/00_config.R")
library(tidyverse)

# Import harmonised eligibility criteria
a_imp_elig <- read_csv("data/elig.csv")

# Get AgD trial set
a_imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Check how many AgD trials are in eligibility (should be all)
a_imp_agd %>% anti_join(a_imp_elig %>% distinct(trial_id))

# Join and isolate trial that is not included
# Select relevant variables
b_elig_agd <- a_imp_agd %>% 
  distinct(trial_id) %>% 
  inner_join(a_imp_elig) %>% 
  select(
    trial_id, 
    contains(c("_age", "bmi", "renal"))
  ) %>% 
  mutate(across(everything(), as.character)) %>% 
  arrange(trial_id) %>% 
  rename(nct_id = trial_id)

# Save
write_csv(b_elig_agd, "../../Attrisk/GitHub/eligibility_criteria.csv")
