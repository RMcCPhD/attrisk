
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

# Collapse variables
# Weighted mean for age, bmi, hba1c and t2dur
# Sum counts for male
b_clp <- b_trttype %>% 
  group_by(trial_id, trttype, arm_label) %>% 
  reframe(
    n = sum(n),
    age = sum(age * n) / sum(n),
    bmi = sum(bmi * n) / sum(n),
    t2dur = sum(t2dur * n) / sum(n),
    hba1c = sum(hba1c * n) / sum(n),
    male = sum(male)
  ) %>% 
  ungroup() %>% 
  distinct()

# Check trials with more than 1 arm that no duplicates remain
b_clp %>%
  count(trial_id) %>%
  filter(n > 1) %>%
  select(trial_id) %>%
  inner_join(b_clp) %>%
  View()

# Join covariates with aggregate groups
c_join <- a_imp_df %>%
  select(trial_id, trttype, arm_n) %>% 
  left_join(
    b_clp %>% select(trial_id, var_n = n, trttype, age:male)
  )

# Get groups with many-to-many joins to resolve these
d_pre <- a_imp_df %>% count(trial_id, trttype, name = "pre")
d_post <- c_join %>% count(trial_id, trttype, name = "post")
d_join <- full_join(d_pre, d_post) %>% mutate(diff = if_else(pre != post, 1, 0))
d_check <- d_join %>% filter(diff == 1) %>% select(trial_id, trttype)

# Manually resolve 4 many-to-many joins (multiple active comparators)
# Easily identified from arm_n
# c_join %>% 
#   inner_join(d_check) %>% 
#   arrange(trial_id) %>% 
#   write_csv("scratch_data/fix_many_joins_vars.csv")

# Resolved using ctgov, arm_n identified between treatment arms
e_fix <- read_csv("scratch_data/fixed_many_joins_vars.csv")

e_fixed <- c_join %>% 
  filter(!(trial_id %in% e_fix$trial_id)) %>% 
  full_join(e_fix)

# Fixed, can now find NAs across variables where possible
# If not, add back in median
e_fixed %>% count(trial_id)
e_fixed %>% reframe(across(everything(), ~ sum(is.na(.))))

write_csv(e_fixed, "processed_data/df_stan2_covars.csv")

