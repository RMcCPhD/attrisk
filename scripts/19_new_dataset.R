
source("scripts/00_config.R")
library(tidyverse)

# Import aggregated trial design features and eligibility
a_imp_df <- read_csv("processed_data/definitive_dataset_with_elig.csv")

# Import imputed and scaled baseline characteristics
a_imp_bl <- read_csv("processed_data/participant_chars_imputed_scaled.csv")

# Join datasets
b_join <- a_imp_bl %>% 
  left_join(
    a_imp_df %>% 
      select(trial_id:source, arm_n_attr, duration:max_bmi_cat, p3:n_regions) %>% 
      group_by(trial_id) %>% 
      slice_head(n = 1) %>% 
      ungroup()
  ) %>% 
  rename(
    arm_n = sample,
    attr = arm_n_attr
  )

# Two trials with incorrect sample sizes in T2D SR repo extraction
# Removing these trials, attrition count bigger than arm_n
c_rm_err <- b_join %>% filter(!(attr > arm_n))

# Final tidying
# Remove hf, some missing
d_tidy <- c_rm_err %>% 
  select(
    source, trial_id:trttype, n, arm_n, attr, duration,
    hf:n_regions, 
    everything()
  ) %>% 
  mutate(across(c(trttype, hf:reg_em), as.integer)) %>% 
  select(-hf)

# Now 344 trials (258 agg, 86 ipd)
d_tidy %>% count(trial_id, source) %>% count(source)

d_tidy %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "na") %>% 
  filter(na != 0)

# Save
# saveRDS(d_tidy, "processed_data/latest_dataset.csv")
