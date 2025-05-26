
# Test for intercept-only model
# Keeping trials where at least 10% attrition happened and duration >=0.5yrs
rm(list = ls())
gc()

library(tidyverse)

a_imp_df <- readRDS("~/Documents/PhD/GitHub/attrisk/processed_data/definitive_dataset.rds")

# Calculate attrition proportion
# Convert duration to weeks
# Filter
b_get <- a_imp_df %>% 
  select(trial_id, arm_n, arm_n_attr, duration) %>% 
  mutate(
    duration = duration * 52,
    attr_pcnt = round(arm_n_attr / arm_n * 100, 1)
  ) %>% 
  filter(attr_pcnt >= 10)

# STAN dataset
g_trial <- list(
  I = nrow(b_get),
  r = b_get$arm_n_attr,
  N = b_get$arm_n,
  dt = b_get$duration
)

saveRDS(g_trial, "mars/intercept_model4/stan_dat.rds")
