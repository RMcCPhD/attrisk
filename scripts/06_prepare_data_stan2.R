
# Preparing STAN data for trial sub-model
# Addition of participant characteristics also in IPD
# Pulled variables that most trials in aact tidied set have (>=80% of trials)
# Age, sex (male), hba1c, bmi, t2dur
# Age, hba1c, bmi and t2dur mean
# Sex (male) as count or percentage
# Also pulled in n for aligning to arms

source("scripts/00_functions_packages_etc.R")

a_imp_df <- readRDS("processed_data/dat_covars.rds")
a_imp_vars <- readRDS("data/ryan_attrisk_vars.rds")

a_imp_df %>% count(trial_id)
a_imp_vars %>% count(trial_id) # 2 less trials, manually find

# Take out grouping variables for joining
b_join_df <- a_imp_df %>% select(trial_id, n, arm_n, trttype)
