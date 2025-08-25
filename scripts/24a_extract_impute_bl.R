
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

# Split by variable
b_split <- split(b_bl_filter, b_bl_filter$variable)
names(b_split)

# Age --------------------------------------------------------------------------

# 334/344 trials
c_age <- b_split$age %>% 
  filter(second_format != "median") %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second)

# Blood pressure  --------------------------------------------------------------

# 89/344 trials
c_sbp <- b_split$blood_pressure %>% 
  filter(category_level == "systolic") %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second)

# BMI --------------------------------------------------------------------------

# 314/344 trials - 311 mean/sd, 3 categorical
c_bmi_df <- b_split$bmi %>% 
  select(trial_id, arm_id, variable, category_level, first_format:second)

# Tidy mean/sd
c_bmi_mean <- c_bmi_df %>% 
  filter(first_format == "mean") %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second)

# Tidy categorical
c_bmi_cat <- c_bmi_df %>% 
  filter(!is.na(category_level)) %>% 
  mutate(
    category = case_match(
      category_level,
      "25 kg/m^2 to 30 kg/m^2" ~ "25-30",
      "30 kg/m^2 and greater" ~ ">=30",
      "< 25 kg/m^2" ~ "<25",
      ">= 25 kg/m^2 - <27 kg/m^2" ~ "25-26.9",
      ">= 27 kg/m^2 - <30 kg/m^2" ~ "27-29.9",
      ">= 30 kg/m^2" ~ ">=30",
      ">=25 kg/m^2" ~ ">=25",
      ">=27 kg/m^2" ~ ">=27",
      ">=30 kg/m^2" ~ ">=30",
      "Less than 25 kg/m^2" ~ "<25"
    ),
    second = second/100
  ) %>% 
  select(trial_id, arm_id, variable, category, first, second) %>% 
  rename(count = first, prop = second)

# T2DM duration ----------------------------------------------------------------

# 275/344 trials
c_t2dur <- b_split$duration_t2dm %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second)

# eGFR -------------------------------------------------------------------------

# 73/344 trials - 72 mean/sd, 1 categorical
c_egfr_df <- b_split$egfr

# Tidy mean/sd
c_egfr_mean <- c_egfr_df %>% 
  filter(first_format == "mean") %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second)

# Tidy categorical
c_egfr_cat <- c_egfr_df %>% 
  filter(!is.na(category_level)) %>% 
  mutate(
    category = case_match(
      category_level,
      "<15 ml/min/1.7m2" ~ "<15",
      ">=15 <30 ml/min/1.7m2" ~ ">=15-<30",
      ">=30 <60 ml/min/1.7m2" ~ ">=30-<60",
      ">=60 <90 ml/min/1.7m2" ~ ">=60-<90",
      ">=90 ml/min/1.7m2" ~ ">90"
    ),
    second = second/100
  ) %>% 
  select(trial_id, arm_id, variable, category, first, second) %>% 
  rename(count = first, prop = second)

# HbA1c ------------------------------------------------------------------------

# 346/344 trials - one missing sd
c_hba1c <- b_split$hba1c %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second)

# HDL cholesterol --------------------------------------------------------------

# 26/344 trials - one has CV, convert this to SD
c_hdl <- b_split$hdl %>% 
  mutate(
    second = case_when(
      second_format == "coeff of variation" ~ first * (second/100),
      TRUE ~ second
    )
  ) %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second)

# Sex --------------------------------------------------------------------------

# 337/344 trials - convert to proportion female
c_sex <- b_split$male %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(count = first, pcnt = second) %>% 
  mutate(
    pcnt = case_when(
      is.na(pcnt) ~ count/169,
      TRUE ~ pcnt/100
    ),
    pcnt = 1 - pcnt,
    variable = "female"
  ) %>% 
  rename(prop = pcnt) %>% 
  select(-count)

# Race -------------------------------------------------------------------------

# 226/344 trials - collapsed categories
c_race <- b_split$race %>% 
  mutate(
    category = case_when(
      grepl("asian", category_level) ~ "asian",
      grepl("black", category_level) ~ "black",
      grepl("white", category_level) ~ "white",
      TRUE ~ "other"
    ),
    second = second/100,
    second = case_when(
      second > 1 ~ 1, 
      TRUE ~ round(second, 2)
    )
  ) %>% 
  select(trial_id, arm_id, variable, category, first, second) %>% 
  rename(count = first, prop = second) %>% 
  select(-count)

# Smoking status ---------------------------------------------------------------

# 8/344 trials
c_smoke <- b_split$smoking_status %>% 
  mutate(
    category = case_match(
      category_level,
      "ex" ~ "former/never",
      "never" ~ "former/never",
      .default = category_level
    ),
    second = second/100
  ) %>% 
  select(trial_id, arm_id, variable, category, first, second) %>% 
  rename(count = first, prop = second)

# Total cholesterol ------------------------------------------------------------

# 26/344 trials - convert coef of variation to sd
c_total <- b_split$total_cholesterol %>% 
  mutate(
    second = case_when(
      second_format == "coeff of variation" ~ first * (second/100),
      TRUE ~ second
    )
  ) %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second)

# Combine ----------------------------------------------------------------------

# 
