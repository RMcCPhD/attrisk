
source("scripts/00_config.R")
library(tidyverse)

# Import aggregated trial design features and eligibility
a_imp_df <- read_csv("processed_data/definitive_dataset_with_elig.csv")

# Import imputed and scaled baseline characteristics
a_imp_bl <- read_csv("processed_data/participant_chars_imputed_scaled.csv")

# Join datasets
b_join <- a_imp_bl %>% 
  left_join(
    a_imp_df %>% 
      select(trial_id:source, arm_n_attr, duration:max_bmi_cat, p3:n_regions) %>% 
      group_by(trial_id) %>% 
      slice_head(n = 1) %>% 
      ungroup()
  ) %>% 
  rename(
    arm_n = sample,
    attr = arm_n_attr
  )

# Two trials with incorrect sample sizes in T2D SR repo extraction
# Removing these trials, attrition count bigger than arm_n
c_rm_err <- b_join %>% filter(!(attr > arm_n))

# Remove hf, some missing
# Remove phase other, none applicable anymore
d_tidy <- c_rm_err %>% 
  select(
    source, trial_id:trttype, n, arm_n, attr, duration,
    hf:n_regions, 
    everything()
  ) %>% 
  mutate(across(c(trttype, hf:reg_em), as.integer))

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
  select(-c(p3, p4, po, bdbl, btrip, bquad, bol, hf, n_regions)) %>% 
  select(source:max_bmi_cat, phase, blinding, everything())

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
g_attr_fixed <- f_attr %>% 
  anti_join(
    read_csv("scratch_data/fixed_many_to_many_na_attr.csv") %>% 
      select(trial_id, arm_id)
  ) %>% 
  full_join(
    read_csv("scratch_data/fixed_many_to_many_na_attr.csv")
  ) %>% 
  arrange(trial_id, trttype)

# Now 344 trials (258 agg, 86 ipd)
g_attr_fixed %>% count(trial_id, source) %>% count(source)

g_attr_fixed %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "na") %>% 
  filter(na != 0)

# Additional covariates and tidying
# Male proportion instead of n
# Add factors for trial size, arm size and number of sites
# Trial size - small (<200), moderate (200-499), large(>=500)
# Arm size - small (<100), moderate (100-250), large (>=250)
# Sites - binary (0 = <100, 1 = >=100)
h_add <- g_attr_fixed %>% 
  mutate(
    male_prop = male_n/arm_n,
    trial_size = case_when(
      n < 200 ~ "small",
      n >= 200 & n < 500 ~ "moderate",
      n >= 500 ~ "large"
    ),
    trial_size = factor(
      trial_size,
      levels = c("small", "moderate", "large")
    ),
    arm_size = case_when(
      arm_n < 100 ~ "small",
      arm_n >= 100 & arm_n < 250 ~ "moderate",
      arm_n >= 250 ~ "large"
    ),
    arm_size = factor(
      arm_size,
      levels = c("small", "moderate", "large")
    ),
    sites = if_else(n_facilities < 100, 0L, 1L)
  ) %>% 
  select(
    source:class, n:attr,
    trttype, duration, trial_size:sites,
    everything()
  ) %>% 
  select(
    source:n_poms, male_prop, age10_mean:egfr_mean
  ) %>% 
  select(-n_facilities) %>% 
  mutate(across(c(trttype, sites:phase, reg_am:reg_em), as.integer))

h_add %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "na") %>% 
  filter(na != 0)

# Save
# saveRDS(h_add, "processed_data/latest_dataset.rds")
