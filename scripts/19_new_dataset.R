
source("scripts/00_config.R")
library(tidyverse)

# Function to quickly check for NAs
check_nas <- function(df) {
  df %>% 
    reframe(across(everything(), ~ sum(is.na(.)))) %>%
    pivot_longer(everything()) %>% 
    filter(value > 0)
}

# Import aggregated trial design features and eligibility
a_imp_df <- read_csv("processed_data/definitive_dataset_with_elig.csv")

# Import imputed and scaled baseline characteristics
a_imp_bl <- read_csv("processed_data/participant_chars_imputed_scaled.csv")

# Import eligibility
a_imp_elig <- read_csv("processed_data/elig_tidy.csv")

# Import trial durations
a_imp_dur <- read_csv("scratch_data/add_duration.csv")

# Prepare eligibility
# Replacing NAs
# BMI NAs = no limits, set to >=20 - <=45
# Age NAs = no upper limit, set to 150
# Add missing trial
b_elig <- a_imp_elig %>% 
  select(trial_id, contains(c("age", "bmi", "egfr"))) %>% 
  mutate(
    age_upper = case_when(
      is.na(age_upper) ~ 150,
      TRUE ~ age_upper
    ),
    bmi_lower = case_when(
      is.na(bmi_lower) ~ "None",
      TRUE ~ bmi_lower
    ),
    bmi_upper = case_when(
      is.na(bmi_upper) ~ "None",
      TRUE ~ bmi_upper
    ),
    egfr_lower = case_when(
      egfr_lower <= 30 ~ "<=30",
      egfr_lower > 30 ~ ">30",
      is.na(egfr_lower) ~ "None"
    ),
    egfr_upper = case_when(
      egfr_upper <= 60 ~ "<=60",
      egfr_upper > 60 ~ ">60",
      is.na(egfr_upper) ~ "None"
    )
  ) %>% 
  bind_rows(
    tibble(
      trial_id = "NCT04170998",
      age_lower = 19,
      age_upper = 150,
      bmi_lower = "18.5",
      bmi_upper = "<40",
      egfr_lower = "None",
      egfr_upper = "None"
    )
  ) %>% 
  distinct()

# Join datasets
b_join <- a_imp_bl %>% 
  left_join(
    a_imp_df %>% 
      select(trial_id:source, arm_n_attr:duration, p3:n_regions) %>% 
      group_by(trial_id) %>% 
      slice_head(n = 1) %>% 
      ungroup()
  ) %>% 
  left_join(b_elig) %>% 
  rename(
    arm_n = sample,
    attr = arm_n_attr
  )

b_join %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything()) %>% 
  filter(value != 0)

# Two trials with incorrect sample sizes in T2D SR repo extraction
# Removing these trials, attrition count bigger than arm_n
b_join %>% filter(attr > arm_n)
c_rm_err <- b_join %>% filter(!(attr > arm_n))

# Arrange variables
d_tidy <- c_rm_err %>% 
  select(
    source, trial_id:trttype, n, arm_n, attr, duration,
    age_lower:egfr_upper, duration:n_poms,
    male_n:egfr_mean
  )

# Factorise binary variables (for models, convert to binary for STAN)
# Remove binary variables for factors
# Remove hf (some missing, not imputing)
e_factors <- d_tidy %>% 
  mutate(
    phase = if_else(p3 == 1L, 0L, 1L),
    blinding = case_when(
      bdbl == 1L ~ "double",
      btrip == 1L ~ "triple",
      bquad == 1L ~ "quadruple",
      bol == 1L ~ "open_label"
    ),
    blinding = as.factor(blinding),
    blinding = fct_relevel(blinding, "open_label")
  ) %>% 
  select(-c(p3, p4, po, bdbl, btrip, bquad, bol)) %>% 
  select(source:egfr_upper, phase:blinding, everything())

# Add attrition counts before they were collapsed
f_attr <- e_factors %>% 
  select(-attr) %>% 
  left_join(
    read_csv("data/trials_n392_attrisk.csv") %>% 
      select(trial_id, atc = atc_new, attr = arm_n_attr) %>% 
      mutate(atc = if_else(atc != "placebo", substr(atc, 1, 5), atc))
  ) %>%
  select(source:arm_n, attr, everything())

# Manually fix where many-to-many were introduced or NA attr
f_attr %>% 
  group_by(trial_id, arm_id) %>% 
  reframe(n = n()) %>% 
  filter(n > 1)

# f_attr %>%
#   group_by(trial_id, arm_id) %>%
#   filter(any(n() > 1) | is.na(attr)) %>%
#   write_csv("scratch_data/fix_many_to_many_na_attr.csv")

# Join fixes
g_fix <- read_csv("scratch_data/fixed_many_to_many_na_attr.csv") %>% 
  select(source:attr) %>% 
  left_join(
    read_csv("scratch_data/fix_many_to_many_na_attr.csv") %>% 
      select(-attr)
  ) %>% 
  distinct()

g_attr_fixed <- f_attr %>% 
  anti_join(
    read_csv("scratch_data/fixed_many_to_many_na_attr.csv") %>% 
      select(trial_id, arm_id)
  ) %>% 
  full_join(g_fix) %>% 
  arrange(trial_id, trttype, attr) %>% 
  group_by(trial_id) %>% 
  fill(egfr_lower:egfr_upper, .direction = "downup") %>% 
  ungroup()

# Now 344 trials (258 agg, 86 ipd)
g_attr_fixed %>% count(trial_id, source) %>% count(source)

g_attr_fixed %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "na") %>% 
  filter(na != 0)

# Final tidying
# Swap male n for male prop
# Format variables (lost by fixing in imported csvs)
# Number of POMs as binary, not much info beyond there being 1 or not
# Binary variables for age limits based on scatter plots in 20_model_shape.R
# Age lower: 0 = 18, 1 = not 18
# Age upper: 0 = 80, 1 = no limit
# Min bmi: 0 = Underweight/healthy (<25), 1 = Obese and above (>=25)
# Max bmi: 0 = Obese (<40), 1 = Morbidly obese (>=40)
# Min egfr: 0 = None, 1 = <=30, 2 = >30
# Max egfr: 0 = None, 1 = <=60, 2 = >60
# Reduced ordinal for number of arms and facilities
# 1 - 2 arms, 2 - 3 arms, 3 - 4 arms, 4 - >4 arms
# Log-transformed number of facilities and sample size

# Update 07.10.25
# Collapsing sample size into categorical variable (i.e. small/med/large)
# Non-sensical otherwise - users would not have a definitive number, dependent
# on number of things including recruitment targets
# Sample: Small (<250), Medium (250 - <500), Large (>=500)
g_attr_fixed %>%
  ggplot(aes(x = n)) +
  geom_histogram() +
  scale_x_continuous(n.breaks = 10) +
  theme_classic()

summary(g_attr_fixed$n)

g_attr_fixed %>%
  ggplot(aes(x = n_facilities)) +
  geom_histogram() +
  scale_x_continuous(n.breaks = 10) +
  theme_classic()

summary(g_attr_fixed$n_facilities)
#
#
# N facilities also collapsed into categorical variable (integers)
# Small (<50), Medium (50-<100), Large (>=100)

h_tidy <- g_attr_fixed %>% 
  rename(
    age_lower_dbl = age_lower,
    age_upper_dbl = age_upper
  ) %>% 
  mutate(
    male_prop = male_n/arm_n,
    across(c(trttype, phase, reg_am:reg_em), as.integer),
    across(c(bmi_lower:bmi_upper, blinding), as.factor),
    blinding = fct_relevel(blinding, "open_label"),
    across(c(n_arms, n_poms), as.factor),
    one_pom = if_else(n_poms == 1, 0L, 1L),
    age_lower = if_else(age_lower_dbl == 18, 0L, 1L),
    age_upper = if_else(age_upper_dbl <= 80, 0L, 1L),
    min_bmi = case_when(
      bmi_lower == "None" ~ 0L,
      as.numeric(str_extract(bmi_lower, "[0-9]+")) <= 25 ~ 1L,
      TRUE ~ 2L
    ),
    max_bmi = case_when(
      bmi_upper == "None" ~ 0L,
      as.numeric(str_extract(bmi_upper, "[0-9]+")) <= 40 ~ 1L,
      TRUE ~ 2L
    ),
    min_egfr = case_when(
      egfr_lower == "None" ~ 0L,
      egfr_lower == "<=30" ~ 1L,
      TRUE ~ 2L
    ),
    max_egfr = case_when(
      egfr_upper == "None" ~ 0L,
      egfr_upper == "<=60" ~ 1L,
      TRUE ~ 2L
    ),
    n_arms = case_when(
      n_arms == 2 ~ 0L,
      n_arms == 3 ~ 1L,
      n_arms == 4 ~ 2L,
      TRUE ~ 3L
    ),
    sample = case_when(
      n < 500 ~ 0L,
      n >= 500 & n < 1000 ~ 1L,
      n >= 1000 ~ 2L
    ),
    n_facilities = case_when(
      n_facilities < 50 ~ 0L,
      n_facilities >= 50 & n_facilities < 100 ~ 1L,
      n_facilities >= 100 ~ 2L
    ),
    # n_facilities = log(n_facilities)
    # n_log = n,
    # across(c(n_log, n_facilities), ~ log(.)),
    # across(c(hba1c_mean, sbp_mean:egfr_mean), ~ log(.))
  ) %>% 
  select(
    source:age_upper_dbl, age_lower:age_upper, bmi_lower:bmi_upper, 
    min_bmi:max_egfr, sample, phase:n_facilities, one_pom
    # male_prop, age10_mean:egfr_mean
  ) %>% 
  select(source:class, n:attr, trttype, everything())

h_tidy %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "na") %>% 
  filter(na != 0)

# Collapse by arm --------------------------------------------------------------

# Isolate counts and baseline characteristics
# Collapse multiarms into one (0 - placebo, 1 - ac, 2 - placebo/ac)
i_vars <- h_tidy %>% 
  select(trial_id, atc:trttype) %>% 
  # select(trial_id, atc:trttype, male_prop:egfr_mean) %>% 
  group_by(trial_id) %>% 
  mutate(
    trttype = case_when(
      length(unique(trttype)) > 1 ~ 2L,
      TRUE ~ trttype
    )
  ) %>% 
  ungroup()

# Collapse arm counts
i_clp_arms <- i_vars %>% 
  select(trial_id:trttype) %>% 
  group_by(trial_id, class) %>% 
  mutate(arm_n = sum(arm_n)) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(
    n = if_else(trial_id == "NCT02897349", 206, n),
    arm_n = if_else(trial_id == "NCT02897349", 102, arm_n),
    attr = if_else(trial_id == "NCT00313313", 82, attr)
  )

# Collapse counts by new trttype
i_clp_trttype <- i_clp_arms %>% 
  select(-c(atc, class)) %>% 
  group_by(trial_id) %>% 
  mutate(
    trttype = case_when(
      length(unique(trttype)) != 1 ~ 2L, 
      TRUE ~ trttype
    )
  ) %>% 
  group_by(trial_id, trttype) %>% 
  mutate(across(arm_n:attr, ~ sum(.))) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(attr_prop = attr/arm_n)

# Collapse baseline characteristics
i_clp_bl <- i_vars %>% 
  group_by(trial_id, trttype) %>% 
  reframe(
    total_n = sum(arm_n)
    # male_prop = weighted.mean(male_prop, arm_n),
    # across(c(age10_mean:egfr_mean), ~ weighted.mean(.x, arm_n))
  ) %>% 
  select(-total_n)

# Check distributions before and after collapse
# i_vars %>%
#   pivot_longer(contains("mean")) %>% 
#   ggplot(aes(x = value)) + 
#   geom_density() +
#   facet_wrap(~ name, scales = "free") +
#   theme_bw()
# 
# i_clp_bl %>%
#   pivot_longer(contains("mean")) %>% 
#   ggplot(aes(x = value)) + 
#   geom_density() +
#   facet_wrap(~ name, scales = "free") +
#   theme_bw()

# Combine collapsed rows
i_join <- i_clp_trttype %>% select(-trttype) %>% left_join(i_clp_bl)
check_nas(i_join)

# Join with trial design features
j_join <- i_join %>% 
  left_join(
    h_tidy %>% 
      select(-c(arm_id, atc, class, n:trttype)) %>%
      # select(-c(arm_id, atc, class, n:trttype, male_prop:egfr_mean)) %>% 
      distinct()
  ) %>% 
  select(
    source, trial_id:attr, attr_prop, trttype, duration:one_pom, everything()
  )

j_join %>% count(trial_id, source) %>% count(source)
check_nas(j_join)

# Add duration
k_dur <- j_join %>% 
  left_join(a_imp_dur %>% rename(dur_new = duration)) %>% 
  mutate(
    duration = dur_new/52,
    duration = case_when(
      is.na(duration) & trial_id == "NCT02730377" ~ 2,
      is.na(duration) & trial_id == "NCT02911948" ~ 0.5,
      TRUE ~ duration
    ),
    n_log = log(n)
  ) %>% 
  select(-dur_new)

k_dur %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "na") %>% 
  filter(na != 0)

# Fix categories for three trials where duplicated introduced after
# updating treatment type
l_fix <- k_dur %>% 
  mutate(
    bmi_lower = case_when(
      trial_id %in% c("NCT00103857", "NCT00798161", "NCT01076088") & bmi_lower == "None" ~ NA,
      TRUE ~ bmi_lower
    ),
    bmi_upper = case_when(
      trial_id %in% c("NCT00103857", "NCT00798161", "NCT01076088") & bmi_upper == "None" ~ NA,
      TRUE ~ bmi_upper
    )
  ) %>% 
  filter(!is.na(bmi_lower) & !is.na(bmi_upper))

# l_fix %>% count(trial_id) %>% filter(n > 1)

# Save
saveRDS(l_fix, "processed_data/latest_dataset.rds")
