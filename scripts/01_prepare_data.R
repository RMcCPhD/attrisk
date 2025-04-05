
# Tidy up dataset and get only the placebo or ac arms
# Ends up with 350 trials, some only had experimental (i.e. comp. dpp4 vs glp1)
source("scripts/00_functions_packages_etc.R")

a_imp_atc <- read_csv("created_metadata/atc.csv")
a_imp_df <- read_csv("data/trials_n392_attrisk.csv")

# Remove glp1, sglt2 and dpp4 arms
# 350 trials (267 agg, 83 ipd)
b_id_arms <- a_imp_df %>% 
  select(-c(unq_drugs, atc, drug_name:dose)) %>% 
  mutate(
    atc_short = case_when(
      atc_new == "placebo" ~ "placebo",
      atc_new != "placebo" ~ substr(atc_new, 1, 5)
    )
  ) %>%
  select(-atc_new) %>% 
  left_join(
    a_imp_atc %>% 
      select(class_short, atc_short = code_short) %>% 
      distinct()
  ) %>% 
  filter(!(class_short %in% c("dpp4", "glp1", "sglt2"))) %>% 
  mutate(class_short = if_else(
    is.na(class_short) & atc_short == "placebo", 
    "placebo", 
    class_short
  ))

# Identify placebo and active comparator arms
b_id_type <- b_id_arms %>% 
  mutate(type = case_when(
    class_short == "placebo" ~ class_short,
    is.na(class_short) ~ "ac",
    TRUE ~ "ac"
  )) %>% 
  select(trial_id:arm_n_attr, type)

# Save (remove ipd when trials are out)
write_csv(b_id_type, "processed_data/df_counts_arms.csv")
