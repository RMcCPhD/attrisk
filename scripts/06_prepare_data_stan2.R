
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

# write_csv(e_fixed, "processed_data/df_stan2_covars.csv")

# Import with manual additions for final tidying
f_imp <- read_csv("processed_data/df_stan2_covars_manually_added.csv")
f_imp %>% count(trial_id)

# Remove arms that aren't active comparators
# Labelled as A10BX in JAMA data but are actually GLP-1s (tirzepatide)
f_rem_exp <- f_imp %>% 
  filter(
    !grepl("this is not an active comparator", note),
    !grepl("remove", note)
  ) %>% 
  select(trial_id:male)

f_rem_exp %>% count(trial_id) # 347 trials after removing tirzepatide
f_rem_exp %>% reframe(across(everything(), ~ sum(is.na(.))))

# Impute missingness -----------------------------------------------------------

library(mice)

# Split dataset by trttype
g_split <- split(f_rem_exp, f_rem_exp$trttype)

# Age (only in placebo)
g_age_vars <- g_split$`1` %>% select(trial_id, arm_n, age) %>% distinct()
g_imp_age <- mice(g_age_vars[, c("arm_n", "age")], m = 5, method = "pmm", seed = 123)
g_imp_age_cmp <- complete(g_imp_age, 1)
g_imp_age_cmp <- g_imp_age_cmp %>% cbind(trial_id = g_age_vars$trial_id)

# BMI
g_bmi0_vars <- g_split$`0` %>% select(trial_id, arm_n, bmi) %>% distinct()
g_imp_bmi0 <- mice(g_bmi0_vars[, c("arm_n", "bmi")], m = 5, method = "pmm", seed = 123)
g_imp_bmi0_cmp <- complete(g_imp_bmi0, 1)
g_imp_bmi0_cmp <- g_imp_bmi0_cmp %>% cbind(trial_id = g_bmi0_vars$trial_id)

g_bmi1_vars <- g_split$`1` %>% select(trial_id, arm_n, bmi) %>% distinct()
g_imp_bmi1 <- mice(g_bmi1_vars[, c("arm_n", "bmi")], m = 5, method = "pmm", seed = 123)
g_imp_bmi1_cmp <- complete(g_imp_bmi1, 1)
g_imp_bmi1_cmp <- g_imp_bmi1_cmp %>% cbind(trial_id = g_bmi1_vars$trial_id)

# T2d duration
g_t2dur0_vars <- g_split$`0` %>% select(trial_id, arm_n, t2dur) %>% distinct()
g_imp_t2dur0 <- mice(g_t2dur0_vars[, c("arm_n", "t2dur")], m = 5, method = "pmm", seed = 123)
g_imp_t2dur0_cmp <- complete(g_imp_t2dur0, 1)
g_imp_t2dur0_cmp <- g_imp_t2dur0_cmp %>% cbind(trial_id = g_t2dur0_vars$trial_id)

g_t2dur1_vars <- g_split$`1` %>% select(trial_id, arm_n, t2dur) %>% distinct()
g_imp_t2dur1 <- mice(g_t2dur1_vars[, c("arm_n", "t2dur")], m = 5, method = "pmm", seed = 123)
g_imp_t2dur1_cmp <- complete(g_imp_t2dur1, 1)
g_imp_t2dur1_cmp <- g_imp_t2dur1_cmp %>% cbind(trial_id = g_t2dur1_vars$trial_id)

# Hba1c
g_hba1c0_vars <- g_split$`0` %>% select(trial_id, arm_n, hba1c) %>% distinct()
g_imp_hba1c0 <- mice(g_hba1c0_vars[, c("arm_n", "hba1c")], m = 5, method = "pmm", seed = 123)
g_imp_hba1c0_cmp <- complete(g_imp_hba1c0, 1)
g_imp_hba1c0_cmp <- g_imp_hba1c0_cmp %>% cbind(trial_id = g_hba1c0_vars$trial_id)

g_hba1c1_vars <- g_split$`1` %>% select(trial_id, arm_n, hba1c) %>% distinct()
g_imp_hba1c1 <- mice(g_hba1c1_vars[, c("arm_n", "hba1c")], m = 5, method = "pmm", seed = 123)
g_imp_hba1c1_cmp <- complete(g_imp_hba1c1, 1)
g_imp_hba1c1_cmp <- g_imp_hba1c1_cmp %>% cbind(trial_id = g_hba1c1_vars$trial_id)

# Male
g_male0_vars <- g_split$`0` %>% select(trial_id, arm_n, male) %>% distinct()
g_imp_male0 <- mice(g_male0_vars[, c("arm_n", "male")], m = 5, method = "pmm", seed = 123)
g_imp_male0_cmp <- complete(g_imp_male0, 1)
g_imp_male0_cmp <- g_imp_male0_cmp %>% cbind(trial_id = g_male0_vars$trial_id)

# Join imputed values
g_imp_join1 <- g_split$`1` %>% 
  select(-c(age)) %>% 
  left_join(g_imp_age_cmp) %>% 
  left_join(
    g_imp_bmi1_cmp, 
    by = c("trial_id", "arm_n"), 
    suffix = c(".old", ".new")
  ) %>% 
  left_join(
    g_imp_t2dur1_cmp, 
    by = c("trial_id", "arm_n"), 
    suffix = c(".old", ".new")
  ) %>% 
  left_join(
    g_imp_hba1c1_cmp, 
    by = c("trial_id", "arm_n"), 
    suffix = c(".old", ".new")
  ) %>% 
  mutate(
    bmi = coalesce(bmi.old, bmi.new),
    t2dur = coalesce(t2dur.old, t2dur.new),
    hba1c = coalesce(hba1c.old, hba1c.new)
  ) %>% 
  select(-contains(c("old", "new")))

g_imp_join0 <- g_split$`0` %>% 
  left_join(
    g_imp_bmi0_cmp, 
    by = c("trial_id", "arm_n"), 
    suffix = c(".old", ".new")
  ) %>% 
  left_join(
    g_imp_t2dur0_cmp, 
    by = c("trial_id", "arm_n"), 
    suffix = c(".old", ".new")
  ) %>% 
  left_join(
    g_imp_hba1c0_cmp, 
    by = c("trial_id", "arm_n"), 
    suffix = c(".old", ".new")
  ) %>% 
  left_join(
    g_imp_male0_cmp, 
    by = c("trial_id", "arm_n"), 
    suffix = c(".old", ".new")
  ) %>% 
  mutate(
    bmi = coalesce(bmi.old, bmi.new),
    t2dur = coalesce(t2dur.old, t2dur.new),
    hba1c = coalesce(hba1c.old, hba1c.new),
    male = coalesce(male.old, male.new)
  ) %>% 
  select(-contains(c("old", "new")))

# Imputation complete, tidying -------------------------------------------------

# Join together
h_join <- g_imp_join1 %>% full_join(g_imp_join0) %>% arrange(trial_id, trttype)
h_join %>% count(trial_id)
h_join %>% reframe(across(everything(), ~ sum(is.na(.))))

# Check counts of each (few changes need manually adjusted in original df)
h_pre <- a_imp_df %>% count(trial_id, trttype, name = "pre")
h_post <- h_join %>% count(trial_id, trttype, name = "post")
h_bind <- full_join(h_pre, h_post) %>% mutate(diff = if_else(pre != post, 1, 0))
h_check <- h_bind %>% filter(diff == 1) %>% select(trial_id, trttype)

# Save a scratch version to remove arms that were collapsed
# i.e. Metformin collapsed into one
# write_csv(a_imp_df %>% filter(trial_id %in% h_join$trial_id), "scratch_data/tmp_trial_df.csv")

# Re-check counts - trttype all match
h_imp <- read_csv("scratch_data/tmp1_trial_df.csv")
h_pre2 <- h_imp %>% count(trial_id, trttype, name = "pre")
h_bind2 <- full_join(h_pre2, h_post) %>% mutate(diff = if_else(pre != post, 1, 0))
h_check2 <- h_bind2 %>% filter(diff == 1) %>% select(trial_id, trttype)

# Find trials with NAs in participant characteristics
i_ids <- h_imp %>% 
  left_join(h_join) %>% 
  group_by(trial_id) %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(-trial_id, names_to = "var", values_to = "na") %>% 
  filter(na != 0) %>% 
  select(trial_id) %>% 
  distinct()

# Filter dataset and manual additions for these trials
# Save to review
# h_imp %>% 
#   filter(trial_id %in% i_ids$trial_id) %>% 
#   arrange(trial_id, trttype) %>% 
#   write_csv("scratch_data/fix_nas_df1.csv")
# 
# f_imp %>% 
#   filter(trial_id %in% i_ids$trial_id) %>% 
#   arrange(trial_id, trttype) %>% 
#   write_csv("scratch_data/fix_nas_df2.csv")

# Import fixes and join
j_fix <- read_csv("scratch_data/fix_nas_df1.csv")

j_join <- h_imp %>% 
  anti_join(j_fix %>% select(trial_id, trttype)) %>% 
  full_join(j_fix) %>% 
  arrange(trial_id, trttype)

j_join %>% count(trial_id) # Same number of trials

# No NAs
j_join %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "na") %>% 
  filter(na != 0)

# Try join again
j_ids <- j_join %>% 
  left_join(h_join) %>% 
  group_by(trial_id) %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(-trial_id, names_to = "var", values_to = "na") %>% 
  filter(na != 0) %>% 
  select(trial_id) %>% 
  distinct()

# Fix arm_n in participant dataset
# h_join %>%
#   filter(trial_id %in% j_ids$trial_id) %>%
#   arrange(trial_id, trttype) %>%
#   write_csv("scratch_data/fix_nas_df3.csv")

# Import and join
k_fix <- read_csv("scratch_data/fix_nas_df3.csv")

k_join <- h_join %>% 
  anti_join(k_fix %>% select(trial_id, trttype)) %>% 
  full_join(k_fix) %>% 
  arrange(trial_id, trttype)

k_join %>% count(trial_id) # Same number of trials

# No NAs (except var_n)
k_join %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "na") %>% 
  filter(na != 0)

# Try join again - no NAs
l_join <- j_join %>% left_join(k_join %>% select(-var_n))
  # group_by(trial_id) %>% 
  # reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  # pivot_longer(-trial_id, names_to = "var", values_to = "na") %>% 
  # filter(na != 0) %>% 
  # select(trial_id) %>% 
  # distinct()

l_join %>% count(trial_id)

# Remove variables no longer needed
# Add proportion male
z_ready <- l_join %>% 
  select(-c(yr:mft_ratio, n_arms:eastern_mediterranean_region)) %>% 
  select(trial_id:arm_n_attr, trttype, duration, male:hba1c, everything()) %>% 
  mutate(male_prop = male/arm_n) %>% 
  select(trial_id:male, male_prop, everything())

# Save
write_csv(z_ready, "processed_data/data_both_models.csv")

# Check low male
z_ready_fixes %>% filter(male_prop <= 0.3324176) %>% 