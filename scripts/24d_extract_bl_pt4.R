
# Pt4 of agg bl extraction
# Checking extraction from 2022 update

source("scripts/00_config.R")
source("scripts/00_paths.R")
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
checkNA(a_imp_set)

# Import pt3 of baseline extraction
a_imp_pt3 <- readRDS("processed_data/agg_bl_third_ext.rds")
checkNA(a_imp_pt3)

# Import read_tabletidier dataset
# No additional data
a_imp_pt4 <- read_csv(z_imp_tt) %>% inner_join(a_imp_pt3 %>% distinct(trial_id))
checkNA(a_imp_pt4)

# Create a lookup of what each trial is missing
a_na_lkp <- a_imp_pt3 %>%
  pivot_longer(age_mu:sbp_d) %>%
  filter(is.na(value)) %>%
  distinct() %>%
  mutate(
    var = case_when(
      grepl("age", name) ~ "age",
      grepl("asian|black|other|white", name) ~ "race",
      grepl("bmi", name) ~ "bmi",
      grepl("current|former", name) ~ "smoke",
      grepl("t2dur_", name) ~ "t2dur",
      grepl("male", name) ~ "sex",
      grepl("sbp", name) ~ "sbp",
      grepl("egfr", name) ~ "egfr",
      grepl("hdl", name) ~ "hdl",
      grepl("total_chol", name) ~ "total_chol",
      grepl("hba1c", name) ~ "hba1c"
    )
  ) %>%
  select(trial_id, var) %>%
  distinct() %>%
  mutate(flag = "1") %>%
  pivot_wider(names_from = var, values_from = flag) %>%
  mutate(
    across(smoke:sex, ~ if_else(is.na(.), "0", .))
  )

# write_csv(a_na_lkp, "scratch_data/bl_missing_pt4.csv")

# Get relevant variables
# Harmonise bl names
b_tt_vars <- a_imp_pt4 %>% 
  select(trial_id, arm_id, variable, category_level, measure_type:second) %>% 
  filter(
    variable %in% c(
      "age", "blood pressure", "blood pressure", "bmi", "duration t2dm",
      "egfr", "hba1c", "hdl", "male", "race", "total cholesterol"
    )
  )

b_tt_vars %>% count(variable)

# Split by variable
b_split <- split(b_tt_vars, b_tt_vars$variable)

# Age --------------------------------------------------------------------------

# 6 trials
c_age_ids <- a_na_lkp %>% filter(age == 1) %>% select(trial_id)

# 0/6
c_age_ext <- b_split$age %>% 
  inner_join(c_age_ids)

# Sex --------------------------------------------------------------------------

# 5 trials
c_sex_ids <- a_na_lkp %>% filter(sex == 1) %>% select(trial_id)

# 0/5
c_sex_ext <- b_split$male %>% 
  inner_join(c_sex_ids)

# Race -------------------------------------------------------------------------

# 308 trials
c_race_ids <- a_na_lkp %>% filter(race == 1) %>% select(trial_id)

# Structured the same as pt1 race data, copying tidying
c_race_ext <- b_split$race %>% 
  inner_join(c_race_ids) %>% 
  mutate(
    second = as.numeric(second),
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
  ungroup()

# 9/132
c_race_ext %>% count(trial_id)

# Create lookup to collapse race levels
# c_race_ext %>%
#   count(category_level) %>%
#   write_csv("scratch_data/race_levels_pt4.csv")

c_race_clp <- c_race_ext %>% 
  left_join(read_csv("scratch_data/race_levels_pt4.csv")) %>% 
  group_by(trial_id, arm_id, race) %>% 
  reframe(prop = sum(second)) %>% 
  group_by(trial_id, arm_id) %>% 
  pivot_wider(names_from = race, values_from = prop) %>% 
  ungroup() %>% 
  mutate(
    total = asian + black + other + white
  )

# Save to manually address rounding errors and check missing against ctgov
# write_csv(c_race_clp, "scratch_data/race_check_pt4.csv")
c_race_fix_imp <- read_csv("scratch_data/race_checked_pt4.csv")

c_race_fix <- c_race_fix_imp %>% 
  select(-note) %>% 
  mutate(total = asian + black + other + white) %>% 
  arrange(trial_id, arm_id)

c_race_new <- a_imp_pt3 %>% 
  select(trial_id, arm_id, asian:white) %>%
  left_join(
    c_race_fix %>% 
      select(-contains("total")) %>% 
      rename_at(vars(black:asian), ~ paste0(., "_new"))
  ) %>% 
  mutate(
    asian = if_else(is.na(asian) & !is.na(asian_new), asian_new, asian),
    black = if_else(is.na(black) & !is.na(black_new), black_new, black),
    other = if_else(is.na(other) & !is.na(other_new), other_new, other),
    white = if_else(is.na(white) & !is.na(white_new), white_new, white)
  ) %>% 
  select(-contains("new"))

# 285 max arms missing race
checkNA(c_race_new)

# BMI --------------------------------------------------------------------------

# 31 trials
c_bmi_ids <- a_na_lkp %>% filter(bmi == 1) %>% select(trial_id)

# 0/31
c_bmi_ext <- b_split$bmi %>% 
  inner_join(c_bmi_ids)

# T2D duration -----------------------------------------------------------------

# 61 trials
c_t2dur_ids <- a_na_lkp %>% filter(t2dur == 1) %>% select(trial_id)

# 1/61 trials
c_t2dur_ext <- b_split$`duration t2dm` %>% 
  inner_join(c_t2dur_ids) %>% 
  separate(second, into = c("1st", "2nd"), sep = ";| ;") %>% 
  mutate(
    across(`1st`:`2nd`, as.numeric),
    mean_t2dur = (`1st` + first + `2nd`) / 3,
    sd_t2dur = (`2nd` - `1st`) / 3
  ) %>% 
  select(trial_id, arm_id, mean_t2dur, sd_t2dur)

c_t2dur_ext %>% count(trial_id)

# Construct t2dur
c_t2dur_new <- a_imp_pt3 %>% 
  select(trial_id, arm_id, t2dur_mu, t2dur_d) %>%
  left_join(c_t2dur_ext %>% mutate(across(contains("t2dur"), as.numeric))) %>% 
  mutate(
    t2dur_mu = if_else(!is.na(t2dur_mu), t2dur_mu, mean_t2dur),
    t2dur_d = if_else(!is.na(t2dur_d), t2dur_d, sd_t2dur)
  ) %>% 
  select(-c(mean_t2dur:sd_t2dur))

# 162 arms missing t2dur (excl. two trials not in aact)
checkNA(c_t2dur_new)

# SBP --------------------------------------------------------------------------

# 255 trials
c_sbp_ids <- a_na_lkp %>% filter(sbp == 1) %>% select(trial_id)

# 0/255 trials
c_sbp_ext <- b_split$`blood pressure` %>% 
  inner_join(c_sbp_ids)

# eGFR -------------------------------------------------------------------------

# 271 trials
c_egfr_ids <- a_na_lkp %>% filter(egfr == 1) %>% select(trial_id)

# 0/271 trials
c_egfr_ext <- b_split$egfr %>% 
  inner_join(c_egfr_ids)

# Total cholesterol ------------------------------------------------------------

# 318 trials
c_total_chol_ids <- a_na_lkp %>% filter(total_chol == 1) %>% select(trial_id)

# 0/318
c_total_chol_ext <- b_split$`total cholesterol` %>% 
  inner_join(c_total_chol_ids)

# HDL cholesterol --------------------------------------------------------------

# 321 trials
c_hdl_ids <- a_na_lkp %>% filter(hdl == 1) %>% select(trial_id)

# 0/318
c_hdl_ext <- b_split$hdl %>% 
  inner_join(c_hdl_ids)

# Smoking status ---------------------------------------------------------------

# 332 trials
c_smoke_ids <- a_na_lkp %>% filter(smoke == 1) %>% select(trial_id)

# 0/322 trials - smoking not in 2022 update set

# Hba1c ------------------------------------------------------------------------

# 35 trials
c_hba1c_ids <- a_na_lkp %>% filter(hba1c == 1) %>% select(trial_id)

c_hba1c_ext <- b_split$hba1c %>% 
  inner_join(c_hba1c_ids) %>% 
  select(trial_id, arm_id, mean_hba1c = first, sd_hba1c = second)

# 13/35
c_hba1c_ext %>% count(trial_id)

# Construct hba1c
c_hba1c_new <- a_imp_pt3 %>% 
  select(trial_id, arm_id, hba1c_mu = mean_hba1c, hba1c_d = sd_hba1c) %>%
  left_join(c_hba1c_ext) %>% 
  mutate(
    across(mean_hba1c:sd_hba1c, as.numeric),
    hba1c_mu = if_else(!is.na(hba1c_mu), hba1c_mu, mean_hba1c),
    hba1c_d = if_else(!is.na(hba1c_d), hba1c_d, sd_hba1c)
  ) %>% 
  select(-c(mean_hba1c:sd_hba1c))

# 43 arms missing hba1c (check those missing sd manually)
checkNA(c_hba1c_new)

# Combine ----------------------------------------------------------------------

# Get dfs
d_get_dfs <- lapply(ls(pattern = "_new"), get)
names(d_get_dfs) <- ls(pattern = "_new")

# Join
d_join <- reduce(d_get_dfs, left_join)

# Join other variables
d_all <- d_join %>% 
  full_join(
    a_imp_pt3 %>%
      select(
        trial_id, arm_id, current, former_never, 
        contains(c("age", "male", "bmi", "sbp", "egfr", "total", "hdl"))
      )
  ) %>% 
  select(
    trial_id:arm_id, contains("age"), male, asian:white, contains("bmi"),
    contains("t2dur"), contains("sbp"), contains("egfr"), contains("total"),
    contains("hdl"), contains("hba1c"), current, former_never
  )


# Check before and after
checkNA(a_imp_pt3) %>% print(n = 24)
checkNA(d_all) %>% print(n = 24)

# Before (AACT, NMA covariates and pre-2022 table tidier)
# Var          before   after
# arm_id            2       2
# age_mu           11      11
# age_d            11      11
# male             10      10
# asian           304     282
# black           307     285
# other           306     284
# white           300     278
# bmi_mu           73      73
# bmi_d            73      73
# mean_hba1c       83      43
# sd_hba1c         90      50
# t2dur_mu        165     162
# t2dur_d         165     162
# sbp_mu          671     671
# sbp_d           671     671
# egfr_mu         711     711
# egfr_d          711     711
# total_chol_mu   849     849
# total_chol_d    849     849
# hdl_mu          852     852
# hdl_d           852     852
# current         884     884
# former_never    884     884

# Save
saveRDS(d_all, "processed_data/agg_bl_fourth_ext.rds")

