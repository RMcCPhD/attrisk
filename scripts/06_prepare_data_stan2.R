
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
a_imp_vars %>% count(trial_id) # 6 less trials, manually find

# Correct trttype in aact vars
b_trttype <- a_imp_vars %>%
  mutate(
    trttype = case_match(
      trttype,
      0 ~ 1L,
      1 ~ 0L
    )
  )

# Join based on grouping variables
a_imp_df %>% filter(trial_id %in% b_trttype$trial_id) %>% count(trial_id, trttype)
b_trttype %>% count(trial_id, trttype)

c_join <- a_imp_df %>% 
  select(trial_id, arm_n, trttype) %>% 
  left_join(
    b_trttype %>% 
      rename(var_arm_n = n) %>% 
      select(trial_id, trttype, arm_label, var_arm_n, age:male)
  )

# Get groups with many-to-many joins to manually resolve these
d_pre <- a_imp_df %>% count(trial_id, trttype, name = "pre")
d_post <- c_join %>% count(trial_id, trttype, name = "post")
d_join <- full_join(d_pre, d_post) %>% mutate(diff = if_else(pre != post, 1, 0))
d_check <- d_join %>% filter(diff == 1) %>% select(trial_id)

# Save to manually address
c_join %>%
  filter(trial_id %in% d_check$trial_id) %>%
  arrange(trial_id) %>% 
  write_csv("scratch_data/fix_many_joins_vars.csv")
