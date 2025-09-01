
# Pt3 of agg bl extraction
# Checking extractions from read_tabletidier

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

# Import first baseline char extraction
# Missingness is two trials not in aact set, manually add at end
a_imp_pt2 <- readRDS("processed_data/agg_bl_second_ext.rds")
checkNA(a_imp_pt2)

# Create a lookup of what each trial is missing
a_na_lkp <- a_imp_pt2 %>%
  pivot_longer(age_mu:total_chol_d) %>%
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
      grepl("sbp", name) ~ "sbp",
      grepl("egfr", name) ~ "egfr",
      grepl("hdl", name) ~ "hdl",
      grepl("chol", name) ~ "total_chol"
    )
  ) %>%
  select(trial_id, var) %>%
  distinct() %>%
  mutate(flag = "1") %>%
  pivot_wider(names_from = var, values_from = flag) %>%
  mutate(
    across(hdl:sex, ~ if_else(is.na(.), "0", .)),
    hba1c = "1"
  )

# write_csv(a_na_lkp, "scratch_data/bl_missing_pt3.csv")

# Import read_tabletidier baseline data
a_imp_tt <- read_csv("../read_tabletidier/data/final_baseline_data.csv")

# Get relevant variables
# Harmonise bl names
b_tt_vars <- a_imp_tt %>% 
  select(trial_id, arm_id, variable, category_level, measure_type:second) %>% 
  filter(
    variable %in% c(
      "age", "blood pressure", "blood_pressure", "bmi", "duration t2dm",
      "duration_t2dm", "egfr", "hba1c", "hdl", "male", "race", "smoking_status",
      "total cholesterol", "total_cholesterol"
    )
  ) %>% 
  mutate(
    variable = case_match(
      variable,
      "blood pressure" ~ "blood_pressure",
      "duration t2dm" ~ "duration_t2dm",
      "total cholesterol" ~ "total_cholesterol",
      .default = variable
    )
  )

b_tt_vars %>% count(variable)

# Split by variable
b_split <- split(b_tt_vars, b_tt_vars$variable)


# Age --------------------------------------------------------------------------

# 10 trials
c_age_ids <- a_na_lkp %>% filter(age == 1) %>% select(trial_id)

# Approximate mean/sd where median, range and/or iqr is given
c_age_ext <- b_split$age %>% 
  inner_join(c_age_ids) %>% 
  separate(second, into = c("1st", "2nd"), sep = "\\-| ") %>% 
  mutate(
    `1st` = gsub("\\(", "", `1st`),
    `2nd` = gsub("\\)", "", `2nd`),
    across(`1st`:`2nd`, as.numeric),
    mean_age = case_when(
      first_format == "mean" ~ first,
      first_format == "median"
      & second_format == "range" ~ (`1st` + 2 * first + `2nd`) / 4,
      first_format == "median"
      & second_format == "iqr" ~ (`1st` + first + `2nd`) / 3
    ),
    sd_age = case_when(
      second_format == "sd" ~ `1st`,
      first_format == "mean"
      & second_format == "range" ~ (`2nd` - `1st`) / 5,
      first_format == "median"
      & second_format == "range" ~ (`2nd` - `1st`) / 3,
      first_format == "median"
      & second_format == "iqr" ~ (`2nd` - `1st`) / 3
    )
  ) %>% 
  select(trial_id, arm_id, mean_age, sd_age)

# 5/10
c_age_ext %>% count(trial_id)

# Construct age
c_age_new <- a_imp_pt2 %>% 
  select(trial_id, arm_id, age_mu, age_d) %>%
  left_join(c_age_ext) %>% 
  mutate(
    age_mu = if_else(!is.na(age_mu), age_mu, mean_age),
    age_d = if_else(!is.na(age_d), age_d, sd_age)
  ) %>% 
  select(-c(mean_age:sd_age))

# 13 arms missing age (excl. two trials not in aact)
checkNA(c_age_new)

# Sex --------------------------------------------------------------------------

# 5 trials
c_sex_ids <- a_na_lkp %>% filter(sex == 1) %>% select(trial_id)

c_sex_ext <- b_split$male %>% 
  inner_join(c_sex_ids) %>% 
  select(trial_id, arm_id, second) %>% 
  mutate(second = round(as.numeric(second)/100, 2))

# 2/10
c_sex_ext %>% count(trial_id)

# Construct sex
c_sex_new <- a_imp_pt2 %>% 
  select(trial_id, arm_id, male) %>%
  left_join(c_sex_ext) %>% 
  mutate(
    male = if_else(!is.na(male), male, second)
  ) %>% 
  select(-second)

# 8 arms missing sex (excl. two trials not in aact)
checkNA(c_sex_new)

# Race -------------------------------------------------------------------------

# 308 trials
c_race_ids <- a_na_lkp %>% filter(race == 1) %>% select(trial_id)

# Structured the same as pt1 race data, copying tidying
c_race_ext <- b_split$race %>% 
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
  ungroup() %>% 
  inner_join(c_race_ids)

# 182/308
c_race_ext %>% count(trial_id)

# Create lookup to collapse race levels
# c_race_ext %>%
#   count(category_level) %>%
#   write_csv("scratch_data/race_levels_pt3.csv")

c_race_clp <- c_race_ext %>% 
  left_join(read_csv("created_metadata/race_levels_pt3_bl.csv")) %>% 
  group_by(trial_id, arm_id, race) %>% 
  reframe(prop = sum(second)) %>% 
  group_by(trial_id, arm_id) %>% 
  mutate(
    prop = if_else(prop == 1.01, 1, prop),
    total = sum(prop)
  ) %>% 
  pivot_wider(names_from = race, values_from = prop)

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
    total_new = round(sum(value_new, na.rm = TRUE), 2),
    value_new = case_when(
      total_new == 0.99 & name == "white" ~ value_new + 0.01,
      total_new == 0.99 & name == "asian" & value_new == 0.99 ~ 1,
      TRUE ~ value_new
    ),
    total_new = round(sum(value_new, na.rm = TRUE), 2)
  ) %>% 
  ungroup() %>% 
  select(-value) %>% 
  pivot_wider(names_from = name, values_from = value_new)

# Isolate trials with prop <1 (10)
# c_race_lt1 <- c_race_fill %>% filter(total_new < 1)
# c_race_lt1 %>% count(trial_id)
# write_csv(c_race_lt1, "scratch_data/race_check_lt1_pt3.csv")

c_race_fix_imp <- read_csv("scratch_data/race_checked_lt1_pt3.csv")

c_race_fix_lt1 <- c_race_fill %>% 
  anti_join(c_race_fix_imp %>% select(trial_id, arm_id)) %>% 
  full_join(c_race_fix_imp %>% select(-note)) %>% 
  group_by(trial_id, arm_id) %>% 
  mutate(
    total_new = asian + black + other + white,
    total_new = if_else(is.na(total_new) & !is.na(asian), 1, total_new)
  ) %>% 
  arrange(trial_id, arm_id)

c_race_new <- a_imp_pt2 %>% 
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
  mutate(
    across(asian:white, ~ round(., 2)),
    white = case_when(
      asian + black + other + white == 0.99 ~ white + 0.01,
      TRUE ~ white
    ),
    other = case_when(
      asian + black + other + white == 1.01 ~ other - 0.01,
      TRUE ~ other
    )
  )

# 321 max arms missing race (excl. two trials not in aact)
checkNA(c_race_new)

# BMI --------------------------------------------------------------------------

# 33 trials
c_bmi_ids <- a_na_lkp %>% filter(bmi == 1) %>% select(trial_id)

c_bmi_ext <- b_split$bmi %>% 
  inner_join(c_bmi_ids) %>% 
  filter(!(first_format == "count" & second_format == "percentage")) %>% 
  separate(second, into = c("1st", "2nd"), sep = ";| ;") %>% 
  mutate(
    across(`1st`:`2nd`, as.numeric),
    mean_bmi = case_when(
      first_format == "mean" ~ first,
      first_format == "median"
      & second_format == "range" ~ (`1st` + 2 * first + `2nd`) / 4,
      first_format == "median"
      & second_format == "iqr" ~ (`1st` + first + `2nd`) / 3
    ),
    sd_bmi = case_when(
      second_format == "sd" ~ `1st`,
      first_format == "mean"
      & second_format == "range" ~ (`2nd` - `1st`) / 5,
      first_format == "median"
      & second_format == "range" ~ (`2nd` - `1st`) / 3,
      first_format == "median"
      & second_format == "iqr" ~ (`2nd` - `1st`) / 3
    )
  ) %>% 
  select(trial_id, arm_id, mean_bmi, sd_bmi)

# 3/33 trials
c_bmi_ext %>% count(trial_id)

c_bmi_new <- a_imp_pt2 %>% 
  select(trial_id, arm_id, bmi_mu, bmi_d) %>%
  left_join(c_bmi_ext) %>% 
  mutate(
    bmi_mu = if_else(!is.na(bmi_mu), bmi_mu, mean_bmi),
    bmi_d = if_else(!is.na(bmi_d), bmi_d, sd_bmi)
  ) %>% 
  select(-c(mean_bmi:sd_bmi))

# T2D duration -----------------------------------------------------------------

# 65 trials
c_t2dur_ids <- a_na_lkp %>% filter(t2dur == 1) %>% select(trial_id)

# 5/65 trials
c_t2dur_ext <- b_split$duration_t2dm %>% 
  inner_join(c_t2dur_ids) %>% 
  separate(second, into = c("1st", "2nd"), sep = ";| ;") %>% 
  mutate(
    across(`1st`:`2nd`, as.numeric),
    mean_t2dur = case_when(
      first_format == "mean" ~ first,
      first_format == "median"
      & second_format == "range" ~ (`1st` + 2 * first + `2nd`) / 4,
      first_format == "median"
      & second_format == "iqr" ~ (`1st` + first + `2nd`) / 3
    ),
    sd_t2dur = case_when(
      second_format == "sd" ~ `1st`,
      first_format == "mean"
      & second_format == "range" ~ (`2nd` - `1st`) / 5,
      first_format == "median"
      & second_format == "range" ~ (`2nd` - `1st`) / 3,
      first_format == "median"
      & second_format == "iqr" ~ (`2nd` - `1st`) / 3
    )
  ) %>% 
  select(trial_id, arm_id, mean_t2dur, sd_t2dur)

c_t2dur_ext %>% count(trial_id)

# Construct t2dur
c_t2dur_new <- a_imp_pt2 %>% 
  select(trial_id, arm_id, t2dur_mu, t2dur_d) %>%
  left_join(c_t2dur_ext %>% mutate(across(contains("t2dur"), as.numeric))) %>% 
  mutate(
    t2dur_mu = if_else(!is.na(t2dur_mu), t2dur_mu, mean_t2dur),
    t2dur_d = if_else(!is.na(t2dur_d), t2dur_d, sd_t2dur)
  ) %>% 
  select(-c(mean_t2dur:sd_t2dur))

# 191 arms missing t2dur (excl. two trials not in aact)
checkNA(c_t2dur_new)

# SBP --------------------------------------------------------------------------

# 255 trials
c_sbp_ids <- a_na_lkp %>% filter(sbp == 1) %>% select(trial_id)

# 0/255 trials
c_sbp_ext <- b_split$blood_pressure %>% 
  inner_join(c_sbp_ids)

# eGFR -------------------------------------------------------------------------

# 271 trials
c_egfr_ids <- a_na_lkp %>% filter(egfr == 1) %>% select(trial_id)

# 0/271 trials
c_egfr_ext <- b_split$egfr %>% 
  inner_join(c_egfr_ids)

# Total cholesterol ------------------------------------------------------------

# 321 trials
c_total_chol_ids <- a_na_lkp %>% filter(total_chol == 1) %>% select(trial_id)

# 1 but median/iqr and values seem incorrect
c_total_chol_ext <- b_split$total_cholesterol %>% 
  inner_join(c_total_chol_ids)

# HDL cholesterol --------------------------------------------------------------

# 321 trials
c_hdl_ids <- a_na_lkp %>% filter(hdl == 1) %>% select(trial_id)

# 1 but median/iqr and values seem incorrect
c_hdl_ext <- b_split$hdl %>% 
  inner_join(c_hdl_ids)

# Smoking status ---------------------------------------------------------------

# 332 trials
c_smoke_ids <- a_na_lkp %>% filter(smoke == 1) %>% select(trial_id)

# 0/322 trials
c_smoke_ext <- b_split$smoking_status %>% 
  inner_join(c_smoke_ids)

# Hba1c ------------------------------------------------------------------------

# 344 trials
c_hba1c_ids <- a_na_lkp %>% filter(hba1c == 1) %>% select(trial_id)

c_hba1c_ext <- b_split$hba1c %>% 
  inner_join(c_hba1c_ids) %>% 
  separate(second, into = c("1st", "2nd"), sep = ";") %>% 
  mutate(
    across(`1st`:`2nd`, as.numeric),
    mean_hba1c = first,
    sd_hba1c = case_when(
      second_format == "sd" ~ `1st`,
      first_format == "mean"
      & second_format == "range" ~ (`2nd` - `1st`) / 5
    )
  ) %>% 
  select(trial_id, arm_id, mean_hba1c, sd_hba1c)

# 314/344
c_hba1c_ext %>% count(trial_id)

# Construct hba1c
c_hba1c_new <- a_imp_pt2 %>% 
  select(trial_id, arm_id) %>%
  left_join(c_hba1c_ext) %>% 
  distinct()

# 90 arms missing hba1c (excl. two trials not in aact)
checkNA(c_hba1c_new)

# Combine ----------------------------------------------------------------------

# Get dfs
d_get_dfs <- lapply(ls(pattern = "_new"), get)
names(d_get_dfs) <- ls(pattern = "_new")

# Join
d_join <- reduce(d_get_dfs, left_join)

# Join other variables
d_all <- a_imp_pt2 %>%
  select(
    trial_id, arm_id, current, former_never, 
    contains(c("total", "hdl", "egfr", "sbp"))
  )

checkNA(d_all)

# Save
saveRDS(d_all, "processed_data/agg_bl_third_ext.rds")

