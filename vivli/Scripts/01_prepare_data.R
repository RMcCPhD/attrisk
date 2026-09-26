
# Join attrition and covariates for 95 trials in AttRisk set

source("Scripts/00_config.R")
library(tidyverse)

# Import data
a_imp_attr <- readRDS("Data/attrition_t2d_n116.rds")

a_imp_t2d <- readRDS(paste0(
  "V:/From 8697/work_folder/model_covariates/",
  "Scratch_data/covariates_with_imputation1e.Rds"
))

a_imp_t2d2 <- readRDS(paste0(
  "V:/From 8697/work_folder/model_covariates/Scratch_data/",
  "data_for_hba1c_main_models_before_standarising_dropmissing_useless.Rds"
))

a_imp_hba1c <- readRDS(paste0(
  "V:/From 8697/work_folder/model_covariates/",
  "Scratch_data/data_for_hba1c_main_models.Rds"
))

a_imp_comorbs <- readRDS("../Extract_data/Comorbs_conmeds/comorbs_all.rds")
a_imp_conmeds <- readRDS("../Extract_data/Comorbs_conmeds/conmeds_all.rds")

# Import lookups
a_imp_ids <- read_csv("Created_metadata/ipd_ids.csv")
a_imp_lkp <- read_csv("Created_metadata/lookup_trials.csv")

a_imp_lkp_arm <- read_csv(paste0(
  "V:/From 8697/work_folder/model_covariates/",
  "Supporting/arm_lvl_info_as_analysed.csv"
))

a_imp_lkp_arm2 <- read_csv(paste0(
  "V:/From 8697/work_folder/model_covariates/",
  "Supporting/arm_lvl_info_as_analysed_out.csv"
))

# Prepare covariates data ------------------------------------------------------

# Combine conmed and T2DM conmed counts
b_cmbn_comorbs <- a_imp_t2d %>% 
  select(nct_id, id) %>% 
  left_join(
    a_imp_conmeds %>% 
      rename(id = usubjid) %>% 
      rename(conmeds = n_cm)
  ) %>% 
  left_join(
    a_imp_comorbs %>% 
      rename(id = usubjid) %>% 
      select(-arm)
  ) %>% 
  mutate(across(c(conmeds, comorbs), ~ if_else(is.na(.), 0, .)))

# Join covariates (incl. hba1c, comorbidity and conmed counts)
b_covs <- a_imp_t2d %>% 
  select(-c(race, race_fix)) %>% 
  left_join(
    a_imp_hba1c %>% 
      select(sponsor, nct_id, data) %>% 
      unnest(data) %>% 
      select(sponsor, nct_id, id, trtcls5, race, hba1c = value_1)
  ) %>% 
  left_join(
    a_imp_t2d2 %>% 
      unnest(data) %>% 
      ungroup() %>% 
      select(nct_id, id, age10, sex)
  ) %>% 
  left_join(b_cmbn_comorbs) %>% 
  select(
    sponsor, nct_id, id, arm, trtcls5,
    age10, sex, bmi, durn, smoke, race, dbp, sbp,
    chol_hdl, chol_ldl, chol_total,
    hba1c, egfr, comorbs, conmeds
  )

# Check value ranges
# Some anomalous values
# Age minimum is 1.4
# BMI minimum is 14.6 but this is possible
# Diastolic BP minimum is 36
# Systolic BP minimum is 69
# HDL minimum is 0.08
# LDL minimum is 0.03
# Total chol minimum is 0.08
b_covs %>% 
  reframe(
    across(
      c(age10, bmi:durn, dbp:egfr),
      list(
        median = ~ median(.x, na.rm = TRUE),
        min = ~ min(.x, na.rm = TRUE),
        max = ~ max(.x, na.rm = TRUE),
        mu = ~ mean(.x, na.rm = TRUE),
        d = ~ sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}-{.fn}"
    )
  ) %>% 
  pivot_longer(everything()) %>% 
  separate(name, into = c("var", "sum"), sep = "-") %>% 
  pivot_wider(names_from = "var", values_from = "value")

# Check how frequent extreme values are
# All very rare (<1%), product of imputation
# Leave as is
b_covs %>% 
  reframe(
    n = n(),
    age = sum(age10 < 1.8, na.rm = TRUE),
    bmi = sum(bmi < 18, na.rm = TRUE),
    dbp = sum(dbp < 60, na.rm = TRUE),
    sbp = sum(sbp < 110, na.rm = TRUE),
    hdl = sum(chol_hdl < 0.3, na.rm = TRUE),
    ldl = sum(chol_ldl < 0.6, na.rm = TRUE),
    total = sum(chol_total < 2, na.rm = TRUE)
  ) %>% 
  pivot_longer(age:total) %>% 
  mutate(prop = value / n)

# Join attrition data and get AttRisk set --------------------------------------

# Join attrition data - 103 trials; 94,432 participants
b_attr <- b_covs %>% 
  left_join(a_imp_attr %>% select(-c(trtgrp, trttype))) %>% 
  select(-c(study)) %>% 
  select(
    sponsor, nct_id, id, arm, trtcls5, 
    days, weeks, event, category, 
    everything()
  )

# Drop 3 trials without attrition extract - 100 trials; 82,381 participants
b_attr_drp <- b_attr %>% 
  group_by(nct_id) %>% 
  filter(!all(is.na(days))) %>%
  ungroup()

# Drop MACE trials - 96 trials; 57,392 participants
b_mace_drp <- b_attr_drp %>% 
  filter(
    !nct_id %in% c(
      "NCT00968708", "NCT01032629", "NCT01131676",
      "NCT01989754", "NCT02065791", "NCT02465515"
    )
  )

# Check attrition missingness
# Spread across 38 trials, 1-201 (most are <5)
b_mace_drp %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

b_mace_drp %>% 
  filter(is.na(days)) %>% 
  count(nct_id, id) %>% 
  group_by(nct_id) %>% 
  reframe(n = sum(n)) %>% 
  print(n = Inf)

# Drop participants with no attrition data - 96 trials; 56,891 participants
b_drp_na <- b_mace_drp %>% filter(!is.na(days))

# Final preparation ------------------------------------------------------------

# Check missingness
# Need to address missing treatment class
# Drop participants missing age or sex
# Covariates dealt with during imputation steps
b_drp_na %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

# b_drp_na %>% 
#   filter(is.na(trtcls5)) %>% 
#   distinct(arm) %>% 
#   write_csv("Scratch_data/address_missing_trtcls5.csv")

# Join added trtcls5 to dataset
c_add_cls <- b_drp_na %>% 
  left_join(
    read_csv("Scratch_data/address_missing_trtcls5.csv") %>% 
      rename(trtcls5_new = trtcls5)
  ) %>% 
  mutate(trtcls5 = coalesce(trtcls5, trtcls5_new)) %>% 
  select(-trtcls5_new)

c_add_cls %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

# Add treatment class name
c_cls_names <- c_add_cls %>% 
  mutate(
    trtcls = case_match(
      trtcls5,
      "A10A" ~ "insulin",
      "A10BA" ~ "biguanide",
      "A10BB" ~ "sulf",
      "A10BF" ~ "agluc",
      "A10BG" ~ "thia",
      "A10BH" ~ "dpp4",
      "A10BJ" ~ "glp1",
      "A10BK" ~ "sglt2",
      "place" ~ "placebo"
    )
  ) %>% 
  select(
    sponsor, nct_id, id, arm, trtcls5, trtcls, everything()
  )

# Drop participants missing age or sex - 96 trials; 56,643 participants
c_drp_na <- c_cls_names %>% filter(!is.na(age10) & !is.na(sex))

# Manually identify non-experimental arms
# c_drp_na %>% 
#   distinct(nct_id, trtcls5, trtcls) %>% 
#   write_csv("Scratch_data/identify_not_exp_arms.csv")

# Keep only placebo or active arms - 91 trials; 20,976 participants
# Drops trials comparing the same drug class (e.g. glp1 vs glp1)
d_drp_exp <- c_drp_na %>% 
  left_join(
    read_csv("Scratch_data/identify_not_exp_arms.csv") %>% 
      filter(not_exp == 1)
  ) %>% 
  filter(!is.na(not_exp)) %>% 
  select(-not_exp)

# Last missingness check
d_drp_exp %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

# Save
saveRDS(d_drp_exp, "Processed_data/01_attr_covariates_n91.rds")
