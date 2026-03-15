
# Summary statistics
source("scripts/00_config.R")
library(tidyverse)

# Import AgD
a_imp_agd <- readRDS("processed_data/latest_dataset.rds")

# AgD: Overall -----------------------------------------------------------------

# Total participants
sum(a_imp_agd$arm_n)

# Follow-up duration
sum_fu <- a_imp_agd %>% 
  mutate(duration = duration * 52) %>% 
  reframe(
    med_fu = median(duration),
    min_fu = min(duration),
    max_fu = max(duration),
    across(med_fu:max_fu, ~ round(.))
  ) %>% 
  mutate(
    value = paste0(med_fu, " [", min_fu, ", ", max_fu, "]")
  )

# Attrition
sum_attr <- a_imp_agd %>% 
  reframe(
    total = sum(arm_n),
    attr = sum(attr),
    attr_p = (attr / total)*100,
    across(total:attr_p, ~ round(., 1))
  ) %>% 
  mutate(
    value = paste0(attr, " [", attr_p, "%]")
  )

sum_attr_risk <- a_imp_agd %>% 
  reframe(
    med_attr = median(attr_prop),
    min_attr = min(attr_prop),
    max_attr = max(attr_prop),
    across(med_attr:max_attr, ~ round(., 2))
  ) %>% 
  mutate(
    value = paste0(med_attr, " [", min_attr, ", ", max_attr, "]")
  )

# Trial-level attrition counts
sum_attr_trial <- a_imp_agd %>% 
  select(trial_id, arm_n, attr) %>% 
  mutate(attr_p = round(attr / arm_n * 100, 2)) %>% 
  mutate(attr = paste0(attr, " (", attr_p, ")")) %>% 
  select(-attr_p)

# Trial characteristics
sum_trial_chars <- a_imp_agd %>% 
  select(phase, blinding, contains("reg"), n_arms, n_facilities, one_pom) %>% 
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(everything()) %>% 
  group_by(name, value) %>% 
  reframe(n = n()) %>% 
  group_by(name) %>% 
  mutate(p = round(n / sum(n) * 100, 2))

# Save
sum_save <- mget(ls(pattern = "sum_"))
lapply(names(sum_save), function(nm) {
  save_dir <- paste0("outputs/summaries/", nm,)
  write_csv(sum_save[[nm]], save_dir)
})

