
# Imputing missing participant characteristics
# Previously imputed age, bmi, t2d duration, hba1c and male proportion
# This was insufficient - only used arm_n as a covariate, and invalid to impute
# male proportion
# Start fresh in this script for the participant characteristics

source("scripts/00_config.R")
library(tidyverse)
library(mice)

# Import latest version of aggregated dataset
a_imp_df <- read_csv("processed_data/definitive_dataset_with_elig.csv")

# Remove participant characteristics
# 347 trials
b_rm_pct <- a_imp_df %>% select(-c(scale_male_prop:hba1c))

# Import latest version of covariates from t2dSR repo
# Get attrisk trials
# 345 trials have data
c_imp_covs <- read_csv("data/base_dsp.csv") %>% 
  inner_join(
    b_rm_pct %>% 
      select(trial_id) %>% 
      distinct()
  )

# Missing trial overall (differs by covariate) - NCT02730377
b_rm_pct %>% anti_join(c_imp_covs %>% select(trial_id) %>% distinct())

# Tidy characteristics ---------------------------------------------------------

# Take relevant variables
d_vars <- c_imp_covs %>% 
  select(trial_id:arm_id, variable, measure_type:second, category_level)

# 2 age missing measure_type, but age always mean where non-NA
d_vars %>% count(variable)
d_vars %>% count(variable, measure_type)

# 2 age where second measure is labelled "median"
# Majority of male prop is a count, some only percentage
d_vars %>% count(variable, measure_type, first_format, second_format)

# Dealing with each covariate individually, then joining together
# Ignoring ethnicity, two trials have it

# Sample size  - 1 trial missing
e_sample <- d_vars %>% 
  filter(variable == "n") %>% 
  select(trial_id:arm_id, sample = first) %>% 
  group_by(trial_id) %>% 
  mutate(n = sum(sample)) %>% 
  ungroup()

# Male - 7 trials missing
# Two missing percentage, keep n to identify
# Calculate n where missing by joining sample sizes
e_male <- d_vars %>% 
  filter(variable == "male") %>% 
  select(trial_id:arm_id, first, second) %>% 
  rename(
    male_n = first,
    male_pcnt = second
  ) %>% 
  mutate(male_n = if_else(!is.na(male_n) & is.na(male_pcnt), male_n, NA)) %>% 
  left_join(e_sample) %>% 
  mutate(
    male_n = case_when(
      is.na(male_n) & !is.na(sample) ~ round(sample * (male_pcnt/100)),
      TRUE ~ male_n
    )
  ) %>% 
  select(-c(sample, male_pcnt))

# Age - 9 trials missing
# Set median for second measure as NA (not usable, rest are sd)
e_age <- d_vars %>% 
  filter(variable == "age") %>% 
  mutate(
    second = case_when(
      second_format == "median" ~ NA,
      TRUE ~ second
    )
  ) %>% 
  rename(
    age_mean = first,
    age_sd = second
  ) %>% 
  select(trial_id:arm_id, age_mean, age_sd)

# Race - 119 trials missing
# Import lookup from t2d sr repo
e_race_lkp <- read_csv("created_metadata/race_smplfy.csv")

# Calculate n using sample and race percentages
# Simplify names and remove collapsed races (name already harmonised by lookup)
e_race <- d_vars %>% 
  filter(variable == "race") %>% 
  select(trial_id:arm_id, category_level, second) %>% 
  rename(
    race = category_level,
    race_pcnt = second
  ) %>% 
  mutate(race_pcnt = if_else(race_pcnt > 100, 100, round(race_pcnt, 2))) %>% 
  left_join(e_sample) %>% 
  mutate(n = case_when(!is.na(sample) ~ round(sample * (race_pcnt/100)))) %>% 
  left_join(e_race_lkp %>% rename(race = category_level)) %>% 
  select(-sample) %>% 
  filter(!grepl("collapsed", race)) %>% 
  rename(race_hrm = smplfy)

# Join together
f_join <- e_sample %>% 
  left_join(e_male) %>% 
  left_join(e_age) %>% 
  left_join(
    e_race %>% 
      select(-c(race, race_pcnt)) %>% 
      mutate(
        race_hrm = case_when(
          race_hrm %in% c("a_ai_h_pi", "hispanic_latino", "other", "unknown") ~ "other",
          TRUE ~ race_hrm
        )
      ) %>% 
      group_by(trial_id, arm_id, race_hrm) %>% 
      reframe(n = sum(n)) %>% 
      pivot_wider(names_from = "race_hrm", values_from = "n")
  )

# Most missingness in race
f_join %>% reframe(across(everything(), ~ sum(is.na(.))))

# Add arms
# Add class for NCT02911948 (based on ctgov results)
# Filter out NA classes (remainder are totals)
# Flag treatment type
# 345 trials
g_arms <- f_join %>% 
  left_join(read_csv("created_metadata/arm_id_trt.csv")) %>% 
  mutate(
    class = case_when(
      trial_id == "NCT02911948" & male_n == 70 ~ "glp1",
      trial_id == "NCT02911948" & male_n == 63 ~ "insulin",
      TRUE ~ class
    ),
    trttype = case_when(
      class == "placebo" ~ 0L,
      !is.na(class) & !(class %in% c("placebo", "dpp4", "glp1", "sglt2")) ~ 1L,
      TRUE ~ NA
    )
  ) %>% 
  filter(!is.na(trttype)) %>% 
  select(trial_id:arm_id, atc:trttype, everything())

g_arms %>% reframe(across(everything(), ~ sum(is.na(.))))

# Add in two missing trials (NCT02730377, NCT03861052)
# NCT03861052 compared two GLP1s, omitting
a_imp_df %>% select(trial_id) %>% distinct() %>% anti_join(g_arms)

# 346 trials
h_add <- g_arms %>% 
  full_join(read_csv("scratch_data/add_n1_trials_chars.csv")) %>% 
  filter(!is.na(trttype)) %>% 
  arrange(trial_id, trttype)

# Additional characteristics ---------------------------------------------------

# Taking additional characteristics from T2D SR repo data
# Relevant variables
z_vars_add <- c("blood_pressure", "duration_t2dm", "egfr", "hba1c", "smoking_status")

# Get variables for 346 trials (335 in additional set)
# Remove disastolic bp
# Only 15 have smoking status, removing
i_imp_add <- read_csv("data/base_dsp_additional.csv") %>% 
  filter(variable %in% z_vars_add) %>% 
  select(trial_id, arm_id, variable:second) %>% 
  inner_join(h_add %>% select(trial_id, arm_id) %>% distinct()) %>% 
  filter(
    !(variable == "blood_pressure" & category_level == "diastolic"),
    variable != "smoking_status",
    !(variable == "egfr" & measure_type == "categorical")
  )

# Pivot wider
j_pivot <- i_imp_add %>% 
  select(trial_id:variable, first, second) %>% 
  mutate(
    variable = case_match(
      variable,
      "blood_pressure" ~ "sbp",
      "duration_t2dm" ~ "t2dur",
      .default = variable
    )
  ) %>% 
  rename(
    mean = first,
    sd = second
  ) %>% 
  pivot_wider(
    names_from = "variable", 
    values_from = c("mean", "sd"),
    names_glue = "{variable}_{.value}",
    names_vary = "slowest"
  )

# Join with other variables
k_join <- h_add %>% left_join(j_pivot)

k_join %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "na") %>% 
  print(n = 21)

# Save before imputation
# write_csv(k_join, "processed_data/participant_chars_before_imputation.csv")

# Address missingness ----------------------------------------------------------

# Imputing male n and means
l_means <- k_join %>% select(trial_id:male_n, contains("mean"))
l_means %>% reframe(across(everything(), ~ sum(is.na(.))))

# Impute
m_impute <- mice(l_means, m = 20, method = "pmm", seed = 123)
m_values <- complete(m_impute, 1) %>% as_tibble()
m_values %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save before scaling
# write_csv(m_values, "processed_data/participant_chars_before_scaling.csv")

# Compare distributions before and after imputation
plot_dist <- m_values %>% 
  pivot_longer(male_n:egfr_mean, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val)) +
  geom_density() +
  geom_density(
    data = l_means %>% 
      pivot_longer(male_n:egfr_mean, names_to = "var", values_to = "val"),
    aes(x = val),
    colour = "forestgreen", 
    alpha = 0.2
  ) +
  facet_wrap(~ var, scales = "free")

# ggsave(
#   "outputs/imputation_before_after.png",
#   plot_dist,
#   width = 6,
#   height = 6,
#   units = "in"
# )

# Scales (standard deviation of means)
m_scales <- m_values %>% 
  pivot_longer(age_mean:egfr_mean, names_to = "var", values_to = "val") %>% 
  group_by(var) %>% 
  reframe(
    mean = mean(val),
    sd = sd(val)
  )

# write_csv(m_scales, "scratch_data/var_scales.csv")

# Scale baseline characteristics after imputation
# Age scaled in decades
# Scaled t2dur, egfr and sbp by sd of means
n_scale <- m_values %>% 
  pivot_longer(
    c(t2dur_mean, sbp_mean:egfr_mean), 
    names_to = "var", 
    values_to = "val"
  ) %>% 
  group_by(var) %>% 
  mutate(
    var = gsub("_mean", "", var),
    val_scaled = round(val / sd(val), 2)
  ) %>% 
  ungroup() %>% 
  select(-val) %>% 
  rename(mean = val_scaled) %>% 
  pivot_wider(
    names_from = "var", 
    values_from = "mean", 
    names_glue = "{var}_{.value}",
    names_vary = "slowest"
  )

# Join scaled variables to rest of variables
o_join <- m_values %>% 
  select(-c(age_mean:egfr_mean)) %>% 
  left_join(n_scale) %>% 
  mutate(age10_mean = age_mean/10) %>% 
  select(-age_mean) %>% 
  select(trial_id:male_n, age10_mean, everything())

# Check distributions after scaling
o_join %>% 
  pivot_longer(male_n:egfr_mean, names_to = "var", values_to = "val") %>% 
  ggplot(aes(x = val)) +
  geom_density() +
  facet_wrap(~ var, scales = "free")

# Last check of missingness
o_join %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save
# write_csv(o_join, "processed_data/participant_chars_imputed_scaled.csv")
