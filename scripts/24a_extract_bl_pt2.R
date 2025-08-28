
source("scripts/00_config.R")
library(tidyverse)

# Helper: Check missingness
checkNA <- function(df) {
  df %>% 
    reframe(across(everything(), ~ sum(is.na(.)))) %>% 
    pivot_longer(everything()) %>% 
    filter(value > 0)
}

# Import trial set
a_imp_set <- readRDS("processed_data/latest_dataset.rds")

# Import first baseline char extraction
# Missingness is two trials not in aact set, manually add at end
a_imp_pt1 <- readRDS("processed_data/agg_bl_first_ext.rds")
checkNA(a_imp_pt1)

# Create a lookup of what each trial is missing
a_na_lkp <- a_imp_pt1 %>%
  pivot_longer(age_m:former_never) %>%
  filter(is.na(value)) %>%
  distinct() %>%
  mutate(
    var = case_when(
      grepl("age", name) ~ "age",
      grepl("asian|black|other|white", name) ~ "race",
      grepl("bmi", name) ~ "bmi",
      grepl("current|former", name) ~ "smoke",
      grepl("dur_", name) ~ "t2dur",
      grepl("male", name) ~ "sex",
      grepl("blood", name) ~ "sbp",
      grepl("egfr", name) ~ "egfr",
      grepl("hdl", name) ~ "hdl",
      grepl("chol", name) ~ "total_chol"
    )
  ) %>%
  select(trial_id, var) %>%
  distinct() %>%
  mutate(flag = "1") %>%
  pivot_wider(names_from = var, values_from = flag) %>%
  mutate(across(sbp:t2dur, ~ if_else(is.na(.), "0", .)))

# write_csv(a_na_lkp, "scratch_data/bl_missing_pt2.csv")

# Pull in harmonised aact datasets ---------------------------------------------

# Get baseline characteristics for trial set
a_imp_bl <- read_csv("data/base_dsp.csv")
a_imp_bl2 <- read_csv("data/base_dsp_additional.csv")

# Combine baseline datasets
a_imp_bl <- a_imp_bl %>% full_join(a_imp_bl2) %>% arrange(trial_id)

a_bl_set <- a_imp_bl %>% inner_join(a_imp_pt1 %>% distinct(trial_id))
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

# Split by variable
b_split <- split(b_bl_filter, b_bl_filter$variable)
names(b_split)

# Age --------------------------------------------------------------------------

# 149 trials
c_age_ids <- a_na_lkp %>% filter(age == 1) %>% select(trial_id)

# 140/149
c_age_ext <- b_split$age %>% 
  filter(second_format != "median") %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(m = first, sd = second) %>% 
  pivot_wider(names_from = variable, values_from = c(m, sd)) %>% 
  inner_join(c_age_ids)

c_age_ext %>% count(trial_id)

# Construct age
c_age_new <- a_imp_pt1 %>% 
  select(trial_id, arm_id, age_m, age_sd) %>%
  left_join(c_age_ext, suffix = c(".old", ".new")) %>% 
  mutate(
    age_mu = if_else(!is.na(age_m), age_m, m_age),
    age_d = if_else(!is.na(age_sd), age_sd, sd_age)
  ) %>% 
  select(-c(age_m:sd_age))

# 19 arms missing age (excl. two trials not in aact)
checkNA(c_age_new)

# Sex --------------------------------------------------------------------------

# 149 trials
c_sex_ids <- a_na_lkp %>% filter(sex == 1) %>% select(trial_id)

c_sex_ext <- b_split$male %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(count = first, pcnt = second) %>% 
  mutate(
    pcnt = case_when(
      is.na(pcnt) ~ count/169,
      TRUE ~ pcnt/100
    )
  ) %>% 
  rename(male = pcnt) %>% 
  select(-count, -variable) %>% 
  inner_join(c_sex_ids)

# 146/149
c_sex_ext %>% count(trial_id)

# Construct male prop
c_sex_new <- a_imp_pt1 %>% 
  select(trial_id, arm_id, male) %>%
  left_join(c_sex_ext %>% rename(male_new = male)) %>% 
  mutate(male = if_else(is.na(male) & !is.na(male_new), male_new, male)) %>% 
  select(-male_new)

# 8 missing male prop (excl. 2 trials not in aact)
checkNA(c_sex_new)

# Race -------------------------------------------------------------------------

# 195 trials
c_race_ids <- a_na_lkp %>% filter(race == 1) %>% select(trial_id)

# Fix proportions >1
# Collapse duplicated groups (found via pivot_wider, lists created)
# Remove "_collapsed" levels
# Add total proportion
# If total >1, take diff from white (see Scratch_data/race_check_mt1.csv)
c_race_ext <- b_split$race %>% 
  mutate(
    second = if_else(second > 100, 100, second),
    second = round(second/100, 2)
  ) %>% 
  select(trial_id:arm_id, category_level, second) %>% 
  group_by(trial_id, arm_id, category_level) %>% 
  reframe(second = sum(second)) %>% 
  filter(!grepl("collapsed", category_level)) %>% 
  group_by(trial_id, arm_id) %>% 
  mutate(
    total = sum(second),
    second = case_when(
      total > 1 & category_level == "white" ~ second - (total - 1),
      TRUE ~ second
    ),
    total = sum(second)
  ) %>% 
  ungroup() %>% 
  inner_join(c_sex_ids)

# Isolate trials where prop >1 (4)
# None left after addressing above
# c_race_mt1 <- c_race_ext %>% filter(total > 1)
# c_race_mt1 %>% count(trial_id)
# write_csv(c_race_mt1, "scratch_data/race_check_mt1.csv")

# 108/195
c_race_ext %>% count(trial_id)

# Create lookup to collapse race levels
# c_race_ext %>% 
#   count(category_level) %>% 
#   write_csv("scratch_data/race_levels.csv")

# Collapse race levels - asian, black, other, white
# All <=1
c_race_clp <- c_race_ext %>% 
  left_join(read_csv("created_metadata/race_levels_pt2_bl.csv")) %>% 
  group_by(trial_id, arm_id, race) %>% 
  reframe(prop = sum(second)) %>% 
  group_by(trial_id, arm_id) %>% 
  mutate(total = sum(prop)) %>% 
  pivot_wider(names_from = race, values_from = prop)

# Get prop for levels where only one is missing (and diff <1)
# Where prop == 1, replace missing with 0
c_race_fill <- c_race_clp %>% 
  pivot_longer(asian:white) %>% 
  group_by(trial_id, arm_id) %>% 
  mutate(
    value_new = case_when(
      is.na(value) & total == 1 ~ 0,
      sum(is.na(value)) == 1 & is.na(value) & total < 1 ~ 1 - total,
      TRUE ~ value
    ),
    value_new = round(value_new, 2),
    total_new = round(sum(value_new, na.rm = TRUE), 2)
  ) %>% 
  ungroup() %>% 
  select(-value) %>% 
  pivot_wider(names_from = name, values_from = value_new)

# Isolate trials with prop <1 (39)
# c_race_lt1 <- c_race_fill %>% filter(total_new < 1)
# c_race_lt1 %>% count(trial_id)
# write_csv(c_race_lt1, "scratch_data/race_check_lt1.csv")

# Fix levels where none missing and prop < 1
c_race_fix_imp <- read_csv("scratch_data/race_checked_lt1.csv")

c_race_fix_lt1 <- c_race_fill %>% 
  anti_join(c_race_fix_imp %>% select(trial_id, arm_id)) %>% 
  full_join(
    c_race_fix_imp %>% 
      select(trial_id, arm_id, total, contains("new")) %>% 
      rename(total_temp = total_new) %>% 
      rename_with(., ~ gsub("_new", "", .)) %>% 
      rename(total_new = total_temp)
  ) %>% 
  group_by(trial_id, arm_id) %>% 
  mutate(total_new = asian + black + other + white) %>% 
  arrange(trial_id, arm_id)

# Construct race
c_race_new <- a_imp_pt1 %>% 
  select(trial_id, arm_id, asian:white) %>%
  left_join(
    c_race_fix_lt1 %>% 
      select(-contains("total")) %>% 
      rename_at(vars(asian:white), ~ paste0(., "_new"))
  ) %>% 
  mutate(
    asian = if_else(is.na(asian) & !is.na(asian_new), asian_new, asian),
    black = if_else(is.na(black) & !is.na(black_new), black_new, black),
    other = if_else(is.na(other) & !is.na(other_new), other_new, other),
    white = if_else(is.na(white) & !is.na(white_new), white_new, white)
  ) %>% 
  select(-contains("new")) %>% 
  mutate(across(asian:white, ~ round(., 2)))

# Repeat tidying for full dataset
c_race_new_tidy <- c_race_new %>% 
  pivot_longer(asian:white) %>% 
  group_by(trial_id, arm_id) %>% 
  mutate(total = round(sum(value, na.rm = TRUE), 2)) %>% 
  ungroup()

# 84 trials with prop >1
# Fix remainder manually
c_race_new_tidy %>% filter(total > 1) %>% count(trial_id)

c_race_fix_mt1 <- c_race_new_tidy %>% 
  filter(total > 1) %>% 
  group_by(trial_id, arm_id) %>% 
  mutate(
    value_new = case_when(
      total == 1.01 & name == "white" ~ value - 0.01,
      value > 1 & value < 2 & name == "white" ~ value - 1,
      value > 2 & name == "white" ~ value - 2,
      TRUE ~ value
    ),
    total_new = sum(value_new, na.rm = TRUE)
  ) %>% 
  ungroup()

# 12 remain after fixes above
# c_race_fix_mt1 %>%
#   filter(total_new > 1) %>%
#   write_csv("scratch_data/race_check_all_mt1.csv")

# Import prop >1 with fixes
c_race_imp_mt1 <- read_csv("scratch_data/race_checked_all_mt1.csv")

c_race_fixed_mt1 <- c_race_fix_mt1 %>% 
  anti_join(c_race_imp_mt1 %>% select(trial_id, arm_id)) %>% 
  full_join(c_race_imp_mt1) %>%
  mutate(
    value_fixed = case_when(
      trial_id %in% c_race_imp_mt1$trial_id ~ value_fix, 
      TRUE ~ value_new
    )
  ) %>% 
  select(trial_id, arm_id, name, value = value_fixed) %>% 
  group_by(trial_id, arm_id) %>% 
  mutate(total = round(sum(value, na.rm = TRUE), 2)) %>% 
  ungroup()

# 178 trials with prop <1 (178)
c_race_new_tidy %>% filter(total < 1) %>% count(trial_id)

c_race_fix_lt1 <- c_race_new_tidy %>% 
  filter(total < 1) %>% 
  group_by(trial_id, arm_id) %>% 
  mutate(
    value_new = case_when(
      sum(is.na(value)) == 1 & is.na(value) ~ 1 - total,
      total == 0.99 & name == "white" ~ value + 0.01,
      TRUE ~ value
    ),
    total_new = sum(value_new, na.rm = TRUE)
  ) %>% 
  ungroup()

# 20 remain (excluding those that are 0)
# c_race_fix_lt1 %>%
#   filter(total_new < 1 & total_new > 0) %>%
#   write_csv("scratch_data/race_check_all_lt1.csv")

# Import prop <1 with fixes
c_race_imp_lt1 <- read_csv("scratch_data/race_checked_all_lt1.csv")

c_race_fixed_lt1 <- c_race_fix_lt1 %>% 
  anti_join(c_race_imp_lt1 %>% select(trial_id, arm_id)) %>% 
  full_join(c_race_imp_lt1) %>%
  mutate(value_fixed = if_else(!is.na(value_fix), value_fix, value_new)) %>% 
  select(trial_id, arm_id, name, value = value_fixed) %>% 
  group_by(trial_id, arm_id) %>% 
  mutate(
    total = round(sum(value, na.rm = TRUE), 2),
    value = if_else(total == 1.01 & name == "white", value - 0.01, value),
    total = round(sum(value, na.rm = TRUE), 2)
  ) %>% 
  ungroup()


# Join fixes - 241 trials
c_race_new_fix <- c_race_fixed_mt1 %>% 
  full_join(c_race_fixed_lt1) %>% 
  select(-total) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  arrange(trial_id)

c_race_new_fix %>% count(trial_id)
