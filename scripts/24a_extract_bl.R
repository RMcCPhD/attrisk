
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

# Helper function - compute pooled SD when summarising group chars
pooledSD <- function(mu, delta, n) {
  # Get weighted mean
  keep_mu <- is.finite(mu) & is.finite(n)
  mu <- mu[keep_mu]
  n_mu <- n[keep_mu]
  if (length(n_mu) == 0) return(NA_real_)
  N <- sum(n_mu)
  if (N <= 1) return(NA_real_)
  Mu <- sum(n_mu * mu) / N
  
  # Get pooled SD
  keep_delta <- is.finite(delta) & is.finite(n)
  within_sd <- sum( (n[keep_delta] - 1) * (delta[keep_delta]^2) )
  bw_sd <- sum( n_mu * (mu - Mu)^2 )
  sqrt( (within_sd + bw_sd) / (N - 1) )
}

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

# Arm size ---------------------------------------------------------------------

c_arm <- b_bl_vars %>% 
  filter(variable == "n") %>% 
  select(trial_id, arm_id, n = first)

# Age --------------------------------------------------------------------------

# 334/344 trials
c_age_mean <- b_split$age %>% 
  filter(second_format != "median") %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second) %>% 
  pivot_wider(names_from = variable, values_from = c(mean, sd))

# Blood pressure  --------------------------------------------------------------

# 89/344 trials
c_sbp_mean <- b_split$blood_pressure %>% 
  filter(category_level == "systolic") %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  mutate(variable = "sbp") %>% 
  rename(mean = first, sd = second) %>% 
  pivot_wider(names_from = variable, values_from = c(mean, sd))

# BMI --------------------------------------------------------------------------

# 314/344 trials - 311 mean/sd, 3 categorical
c_bmi_df <- b_split$bmi %>% 
  select(trial_id, arm_id, variable, category_level, first_format:second)

# Tidy mean/sd
c_bmi_mean <- c_bmi_df %>% 
  filter(first_format == "mean") %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second) %>% 
  pivot_wider(names_from = variable, values_from = c(mean, sd))

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
  rename(count = first, prop = second) %>% 
  arrange(trial_id, arm_id)

# T2DM duration ----------------------------------------------------------------

# 275/344 trials
c_t2dur_mean <- b_split$duration_t2dm %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  mutate(variable = "t2dur") %>% 
  rename(mean = first, sd = second) %>% 
  pivot_wider(names_from = variable, values_from = c(mean, sd))

# eGFR -------------------------------------------------------------------------

# 73/344 trials - 72 mean/sd, 1 categorical
c_egfr_df <- b_split$egfr

# Tidy mean/sd
c_egfr_mean <- c_egfr_df %>% 
  filter(first_format == "mean") %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second) %>% 
  pivot_wider(names_from = variable, values_from = c(mean, sd))

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
c_hba1c_mean <- b_split$hba1c %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second) %>% 
  pivot_wider(names_from = variable, values_from = c(mean, sd))

# HDL cholesterol --------------------------------------------------------------

# 26/344 trials - one has CV, convert this to SD
c_hdl_mean <- b_split$hdl %>% 
  mutate(
    second = case_when(
      second_format == "coeff of variation" ~ first * (second/100),
      TRUE ~ second
    ),
    variable = "hdl"
  ) %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second) %>% 
  pivot_wider(names_from = variable, values_from = c(mean, sd))

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
  rename(female_prop = pcnt) %>% 
  select(-count, -variable)

# Race -------------------------------------------------------------------------

# 226/344 trials - collapsed categories
# Removed collapsed levels
# Fixing totals >1 caused by proportions
c_race <- b_split$race %>% 
  filter(!grepl("collaps", category_level)) %>% 
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
  select(-variable, -count) %>% 
  group_by(trial_id, arm_id, category) %>% 
  reframe(prop = sum(prop)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = category, values_from = prop) %>% 
  mutate(
    across(asian:white, ~ if_else(. > 1, . - 1, .)),
    across(asian:white, ~ if_else(. > 1, 1, .)),
    across(asian:white, ~ if_else(is.na(.), 0, .)),
    race_total = asian + black + other + white,
    other = case_when(
      race_total == 1.29 ~ other - 0.29,
      race_total == 1.04 ~ other - 0.04,
      race_total == 1.03 ~ other - 0.03,
      race_total == 1.02 ~ other - 0.02,
      race_total == 1.01 ~ other - 0.01,
      race_total == 0.99 ~ other + 0.01,
      TRUE ~ other
    ),
    race_total = asian + black + other + white,
    across(asian:race_total, ~ round(., 2))
  )

# Adding missing to other where prop < 1
c_race_fix <- c_race %>% 
  mutate(
    asian = if_else(race_total == 0.01, 1, asian),
    race_total = asian + black + other + white,
    other = if_else(race_total < 1, other + (1 - race_total), other),
    race_total = asian + black + other + white,
    across(asian:race_total, ~ round(., 2)),
    other = if_else(race_total == 1.01, other - 0.01, other),
    race_total = asian + black + other + white
  )

# Smoking status ---------------------------------------------------------------

# 8/344 trials
c_smoke <- b_split$smoking_status %>% 
  mutate(
    category = case_match(
      category_level,
      "ex" ~ "former_never",
      "never" ~ "former_never",
      .default = category_level
    ),
    second = second/100
  ) %>% 
  select(trial_id, arm_id, variable, category, first, second) %>% 
  rename(count = first, prop = second) %>% 
  select(-variable, -count) %>% 
  group_by(trial_id, arm_id, category) %>% 
  reframe(prop = sum(prop)) %>% 
  pivot_wider(names_from = category, values_from = prop) %>% 
  mutate(
    former_never = if_else(is.na(former_never), 1 - current, former_never),
    across(current:former_never, ~ round(., 2)),
    smoke_total = current + former_never
  )

# Total cholesterol ------------------------------------------------------------

# 26/344 trials - convert coef of variation to sd
c_total <- b_split$total_cholesterol %>% 
  mutate(
    second = case_when(
      second_format == "coeff of variation" ~ first * (second/100),
      TRUE ~ second
    ),
    variable = "total_chol"
  ) %>% 
  select(trial_id, arm_id, variable, first, second) %>% 
  rename(mean = first, sd = second) %>% 
  pivot_wider(names_from = variable, values_from = c(mean, sd))

# Combine ----------------------------------------------------------------------

# List of tidied datasets
d_lst <- list(
  c_age_mean, c_sbp_mean, c_bmi_mean,c_t2dur_mean, c_egfr_mean,
  c_hba1c_mean, c_hdl_mean, c_sex, c_race, c_smoke, c_total
)

# Join datasets together
d_join <- reduce(d_lst, full_join)

# Join with IDs from trial set
d_join_ids <- a_imp_set %>% 
  distinct(trial_id) %>% 
  left_join(d_join) %>% 
  select(-race_total, -smoke_total)

# One trial missing all characteristics (manually add in)
# Rest have varied missingness
# Mostly smoking, cholesterol (both), egfr,
d_join_ids %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  mutate(prop = value/1062) %>% 
  arrange(desc(prop)) %>% 
  print(n = 25)

# Join arms from trial set
d_imp_arms <- read_csv("created_metadata/arm_id_trt.csv")

d_join_arms <- d_join_ids %>% 
  left_join(d_imp_arms) %>% 
  left_join(c_arm) %>% 
  select(trial_id, n, atc:class, everything()) %>% 
  filter(!class %in% c("dpp4", "glp1", "sglt2"), !is.na(class)) %>% 
  select(-arm_id, -atc) %>% 
  mutate(across(mean_age:sd_total_chol, ~ round(., 2))) %>% 
  distinct()

# Assign trttype to match trial set
e_trttype <- d_join_arms %>% 
  select(trial_id, class, n) %>% 
  group_by(trial_id) %>% 
  mutate(
    trttype = case_when(
      length(unique(class)) == 1 & class == "placebo" ~ 0L,
      length(unique(class)) == 1 & class != "placebo" ~ 1L,
      TRUE ~ 2L
    )
  ) %>% 
  ungroup()

# Join trttype to baseline set
e_join_trttype <- d_join_arms %>% 
  cbind(e_trttype %>% select(trttype)) %>% 
  as_tibble() %>% 
  select(trial_id, trttype, everything())

# Collapse variables by trttype
# 342 trials
f_clp_chars <- e_join_trttype %>% 
  group_by(trial_id, trttype) %>% 
  reframe(
    across(contains("mean"), ~ weighted.mean(., n)),
    psd_age = pooledSD(mean_age, sd_age, n),
    psd_sbp = pooledSD(mean_sbp, sd_sbp, n),
    psd_bmi = pooledSD(mean_bmi, sd_bmi, n),
    psd_t2dur = pooledSD(mean_t2dur, sd_t2dur, n),
    psd_egfr = pooledSD(mean_egfr, sd_egfr, n),
    psd_hba1c = pooledSD(mean_hba1c, sd_hba1c, n),
    psd_hdl = pooledSD(mean_hdl, sd_hdl, n),
    psd_total_chol = pooledSD(mean_total_chol, sd_total_chol, n),
    across(female_prop:former_never, ~ weighted.mean(., n)),
  ) %>%
  mutate(across(mean_age:former_never, ~ round(., 2))) %>% 
  select(
    trial_id:trttype, contains("age"), contains("sbp"), contains("bmi"),
    contains("t2dur"), contains("egfr"), contains("hba1c"), contains("hdl"),
    contains("total_chol"), everything()
  ) %>% 
  rename_with(., ~ gsub("psd", "sd", .))

# Re-add missing trials
# Align trttype to trial set for two trials
f_all <- f_clp_chars %>% 
  full_join(
    a_imp_set %>% 
      filter(!trial_id %in% f_clp_chars$trial_id) %>% 
      distinct(trial_id, trttype)
  ) %>% 
  arrange(trial_id) %>% 
  mutate(
    trttype = case_when(
      trial_id %in% c("NCT00676338", "NCT01147627") ~ 1L,
      TRUE ~ trttype
    )
  )

# Check trttype breakdown - matches across trials
f_all %>% count(trial_id, trttype) %>% count(trttype)
a_imp_set %>% count(trial_id, trttype) %>% count(trttype)

# Save
saveRDS(f_all, "processed_data/bl_chars_noimp.rds")
