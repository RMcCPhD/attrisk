
source("Scripts/00_config.R")
library(tidyverse)
library(survival)
library(broom)
library(glmnet)

# Import dataset - 95 trials
a_imp_df <- readRDS("Processed_data/n95_conmeds.rds")

# Get placebo/ac arms - 86 trials
a_df <- a_imp_df %>% 
  select(-trttype) %>% 
  group_by(nct_id) %>% 
  filter(!(trt_class %in% c("glp1", "dpp4", "sglt2"))) %>% 
  ungroup() %>% 
  mutate(
    trttype = if_else(trt_class == "placebo", "plc", "ac"),
    event = as.integer(event)
  ) %>% 
  select(sponsor:trt_class, trttype, everything())

a_df %>% count(nct_id)
a_df %>% count(trt_class, trttype)

# Add person-time
a_df_pt <- a_df %>% 
  group_by(nct_id) %>% 
  mutate(
    pyr = weeks/52,
    mft = sum(pyr)/n()
  ) %>% 
  select(sponsor:weeks, pyr, mft, everything())

# Check events per arm type in trials
a_events <- a_df %>% 
  group_by(nct_id, trttype) %>% 
  reframe(events = sum(event)) %>% 
  mutate(low = if_else(events < 10, 1, 0))

# Remove arms in trials with less than 10 events (test to see if it fits)
# 77 trials
b_filter <- a_events %>% filter(low == 1) %>% select(nct_id, trttype)
b_prep_df <- a_df %>% anti_join(b_filter)
b_prep_df %>% count(nct_id)

# Remove participants missing age and sex
b_age_sex <- b_prep_df %>% filter(!is.na(age10) & !is.na(sex))

# Add race into monoethnic trials
# Redundant - race will not contribute to explaining variability in attrition
# in cases where all participants were in one group
# b_add_race <- b_age_sex %>% 
#   mutate(
#     race = case_when(
#       nct_id %in% c("NCT00654381", "NCT01214239", "NCT01215097", "NCT02453555") ~ "asian",
#       nct_id %in% c("NCT01194830", "NCT02182830") ~ "black_aa",
#       TRUE ~ race
#     )
#   )

# Save
saveRDS(b_age_sex, "Processed_data/n77_cox.rds")
