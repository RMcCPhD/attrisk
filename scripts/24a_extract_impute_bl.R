
# Predictors
# Trttype, smoke, sex, sbp
# Race (asian, black/aa, white, other)
# Number of con meds and comorbidities
# HbA1c, eGFR, duration of T2DM
# Total and HDL cholesterol
# BMI, age10

# Extract and impute first
# Collapse into weighted marginals after
source("scripts/00_config.R")
library(tidyverse)
library(multinma)

# Import trial dataset and aggregated baseline characteristics
a_imp_set <- readRDS("processed_data/latest_dataset.rds")
a_imp_bl <- read_csv("data/base_dsp.csv")
a_imp_bl2 <- read_csv("data/base_dsp_additional.csv")

# Combine baseline datasets
a_imp_bl <- a_imp_bl %>% full_join(a_imp_bl2) %>% arrange(trial_id)

# Get baseline characteristics for trial set
a_bl_set <- a_imp_bl %>% inner_join(a_imp_set %>% select(trial_id))
a_bl_set %>% count(trial_id) # 1 missing

# Select relevant variables
b_bl_vars <- a_bl_set %>% select(trial_id:arm_id, variable:second)
b_bl_vars %>% count(trial_id, variable) %>% count(variable) %>% arrange(n)

# Filter for baseline characteristics among pooled Cox estimates
# Got everything except con meds and comorbidities
b_bl_filter <- b_bl_vars %>% 
  filter(
    !variable %in% c(
      "ethnicity", "heart_failure", "total_hdl_ratio",
      "cardiovascular_disease", "n"
    )
  )
