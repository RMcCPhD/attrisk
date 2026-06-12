
rm(list = ls())
gc()

library(tidyverse)

a_imp_df <- read_csv("processed_data/definitive_dataset.csv")

a_imp_df %>% 
  group_by(trttype) %>% 
  reframe(
    trials = length(unique(trial_id)),
    n = sum(n)
  )

# Total
a_imp_df %>% 
  group_by(trial_id) %>% 
  reframe(n = sum(arm_n)) %>% 
  reframe(
    total = sum(n),
    median = median(n),
    min = min(n),
    max = max(n)
  )

a_imp_df %>% 
  group_by(trial_id) %>% 
  reframe(events = sum(arm_n_attr)) %>% 
  reframe(
    total = sum(events),
    median = median(events),
    min = min(events),
    max = max(events)
  )

a_imp_df %>% 
  reframe(
    median = median(duration) * 52,
    min = min(duration) * 52,
    max = max(duration) * 52
  )

# Arms
a_imp_df %>% 
  group_by(trttype, trial_id) %>% 
  reframe(n = sum(arm_n)) %>% 
  group_by(trttype) %>% 
  reframe(
    total = sum(n),
    median = median(n),
    min = min(n),
    max = max(n)
  )

a_imp_df %>% 
  group_by(trttype, trial_id) %>% 
  reframe(events = sum(arm_n_attr)) %>% 
  group_by(trttype) %>% 
  reframe(
    total = sum(events),
    median = median(events),
    min = min(events),
    max = max(events)
  )

a_imp_df %>% 
  group_by(trttype) %>% 
  reframe(
    median = median(duration) * 52,
    min = min(duration) * 52,
    max = max(duration) * 52
  )

