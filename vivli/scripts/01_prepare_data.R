
source("Scripts/00_config.R")
library(tidyverse)

# Data
a_imp_attr <- readRDS("Data/attrition_t2d_n116.rds")
a_imp_t2d <- readRDS("V:/From 8697/work_folder/model_covariates/Scratch_data/data_for_hba1c_main_models.Rds")

# Lookups
a_imp_ids <- read_csv("Created_metadata/ipd_ids.csv")
a_imp_lkp <- read_csv("Created_metadata/lookup_trials.csv")
a_imp_lkp_arm <- read_csv("V:/From 8697/work_folder/model_covariates/Supporting/arm_lvl_info_as_analysed.csv")
a_imp_lkp_arm2 <- read_csv("V:/From 8697/work_folder/model_covariates/Supporting/arm_lvl_info_as_analysed_out.csv")

# Get 95 trials in attrition set
b_n95_attr <- a_imp_attr %>% 
  inner_join(a_imp_ids %>% rename(nct_id = trial_id)) %>% 
  mutate(
    trttype = case_match(
      trttype,
      "Plc/AC" ~ "Plc",
      "Placebo" ~ "Plc",
      .default = trttype
    )
  )

# Get t2d data for 95 trials
b_n95_t2d <- a_imp_t2d %>% 
  select(sponsor:nct_id, data) %>% 
  unnest(data) %>% 
  filter(nct_id %in% b_n95_attr$nct_id) %>% 
  select(
    sponsor:arm_id_unq, trtcls4:trtcls5, age10:tot_hdl_ratio_log, cursmoke, race
  )

# Check NAs - some in trtcls4 and trtcls5
b_n95_t2d %>%
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything(), names_to = "var", values_to = "na") %>% 
  filter(na != 0)

# NAs have ATC coding, can add these in
b_n95_t2d %>% filter(is.na(trtcls5)) %>% count(arm_lvl)

# Add in missing trtcls5
b_add_cls <- b_n95_t2d %>% 
  mutate(
    trtcls4 = case_when(is.na(trtcls4) ~ substr(arm_lvl, 1, 4), TRUE ~ trtcls4),
    trtcls5 = case_when(is.na(trtcls5) ~ substr(arm_lvl, 1, 5), TRUE ~ trtcls5)
  )

b_add_cls %>% filter(is.na(trtcls5)) %>% count(arm_lvl)

# Join attrition and covariates
c_attr_t2d <- b_n95_attr %>% inner_join(b_add_cls)
c_attr_t2d %>% reframe(across(days:category, ~ sum(is.na(.))))

# Add treatment classes
d_join_classes <- c_attr_t2d %>% 
  mutate(
    trt_class = case_match(
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
  )

# Order variables and last missingness check
d_join_classes %>% reframe(across(everything(), ~ sum(is.na(.))))

e_order <- d_join_classes %>% 
  select(
    sponsor:trttype, trt_class, arm_lvl, days:category,
    age10:sex, race, cursmoke, durn:tot_hdl_ratio_log
  ) %>% 
  select(-arm_f) %>% 
  arrange(sponsor, nct_id, days)

# Save
saveRDS(e_order, "Processed_data/n95.rds")
