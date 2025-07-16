
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

# Prepare eligibility
# Replacing NAs
# BMI NAs = no limits, set to >=20 - <=45
# Age NAs = no upper limit, set to 150
# Add missing trial
b_elig <- a_imp_elig %>% 
  select(trial_id, contains(c("age", "bmi"))) %>% 
  mutate(
    age_upper = case_when(
      is.na(age_upper) ~ 150,
      TRUE ~ age_upper
    ),
    bmi_lower = case_when(
      is.na(bmi_lower) ~ ">=20",
      TRUE ~ bmi_lower
    ),
    bmi_upper = case_when(
      is.na(bmi_upper) ~ "<=45",
      TRUE ~ bmi_upper
    )
  ) %>% 
  bind_rows(
    tibble(
      trial_id = "NCT04170998",
      age_lower = 19,
      age_upper = 150,
      bmi_lower = "18.5",
      bmi_upper = "<40"
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
    age_lower:bmi_upper, duration:n_poms,
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
  select(source:bmi_upper, phase:blinding, everything())

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
  arrange(trial_id, trttype, attr)

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
# Reduced ordinal for number of arms and facilities
# 1 - 2 arms, 2 - 3 arms, 3 - 4 arms, 4 - >4 arms
# Log-transformed number of facilities
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
    age_upper = if_else(age_upper_dbl == 80, 0L, 1L),
    min_bmi = case_when(
      as.numeric(str_extract(bmi_lower, "[0-9]+")) <= 25 ~ 0L,
      TRUE ~ 1L
    ),
    max_bmi = case_when(
      as.numeric(str_extract(bmi_upper, "[0-9]+")) <= 40 ~ 0L,
      TRUE ~ 1L
    ),
    n_arms = case_when(
      n_arms == 2 ~ 1L,
      n_arms == 3 ~ 2L,
      n_arms == 4 ~ 3L,
      TRUE ~ 4L
    ),
    n_facilities = log(n_facilities),
    across(c(hba1c_mean, sbp_mean:egfr_mean), ~ log(.))
  ) %>% 
  select(
    source:age_upper_dbl, age_lower:age_upper, bmi_lower:bmi_upper, 
    min_bmi:max_bmi, phase:n_facilities, one_pom,
    male_prop, age10_mean:egfr_mean
  ) %>% 
  select(source:class, n:attr, trttype, everything())

h_tidy %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "na") %>% 
  filter(na != 0)

# Collapse by arm --------------------------------------------------------------

# Isolate counts and baseline characteristics
i_vars <- h_tidy %>% select(trial_id, atc:trttype, male_prop:egfr_mean)

# Collapse arm counts
i_clp_arms <- i_vars %>% 
  select(trial_id:trttype) %>% 
  group_by(trial_id, class) %>% 
  mutate(arm_n = sum(arm_n)) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(
    attr = if_else(trial_id == "NCT00313313", 82, attr),
    attr_prop = attr/arm_n
  )

# Collapse baseline characteristics
i_clp_bl <- i_vars %>% 
  group_by(trial_id, class) %>% 
  reframe(
    total_n = sum(arm_n),
    male_prop = weighted.mean(male_prop, arm_n),
    across(c(age10_mean:egfr_mean), ~ weighted.mean(.x, arm_n))
  ) %>% 
  select(-total_n)

# Check distributions before and after collapse
i_vars %>%
  pivot_longer(contains("mean")) %>% 
  ggplot(aes(x = value)) + 
  geom_density() +
  facet_wrap(~ name, scales = "free") +
  theme_bw()

i_clp_bl %>%
  pivot_longer(contains("mean")) %>% 
  ggplot(aes(x = value)) + 
  geom_density() +
  facet_wrap(~ name, scales = "free") +
  theme_bw()

# Combine collapsed rows
i_join <- i_clp_arms %>% left_join(i_clp_bl)
check_nas(i_join)

# Join with trial design features
j_join <- h_tidy %>% 
  select(source:n, trttype:one_pom) %>% 
  select(-arm_id) %>% 
  distinct() %>% 
  left_join(i_join) %>% 
  select(source:n, arm_n:attr_prop, trttype:one_pom, everything())

# Save
saveRDS(j_join, "processed_data/latest_dataset.rds")
