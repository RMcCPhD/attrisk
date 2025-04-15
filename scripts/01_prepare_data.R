
# Tidy up dataset and get only the placebo or ac arms
# Ends up with 350 trials, some only had experimental (i.e. comp. dpp4 vs glp1)
source("scripts/00_functions_packages_etc.R")

a_imp_atc <- read_csv("created_metadata/atc.csv")
a_imp_df <- read_csv("data/trials_n392_attrisk.csv")
a_imp_ipd <- read_csv("data/attr_n86.csv")
a_imp_fu <- read_csv("data/fu_ratio_n86.csv")

# Remove glp1, sglt2 and dpp4 arms
# Add treatment type
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
  mutate(
    class_short = if_else(
      is.na(class_short) & atc_short == "placebo", 
      "placebo", 
      class_short
    ),
    trttype = case_when(
      class_short == "placebo" ~ "plc",
      is.na(class_short) ~ "ac",
      TRUE ~ "ac"
    ))

# 350 trials (267 agg, 83 ipd)
b_id_arms %>% count(trial_id, source) %>% count(source)

# Prepare IPD counts
c_prp_ipd <- a_imp_ipd %>% 
  mutate(collapsed = strsplit(collapsed, " AND ")) %>% 
  unnest(collapsed) %>% 
  select(-length) %>% 
  separate(
    collapsed,
    into = c("trttype", "n", "arm_n_attr", "compl"),
    sep = " "
  ) %>% 
  mutate(across(trttype:compl, ~ gsub("\\b\\w+=", "", .))) %>% 
  mutate(
    source = "ipd",
    across(n:arm_n_attr, as.numeric)
  ) %>% 
  rename(
    trial_id = nct_id,
    arm_n = n
  )

# 86 trials
c_prp_ipd %>% count(trial_id)

# All there (83/86)
sum(unique(c_prp_ipd$trial_id) %in% b_id_arms$trial_id[b_id_arms$source == "ipd"])

c_add_ipd <- b_id_arms %>% 
  filter(!trial_id %in% c_prp_ipd$trial_id) %>% 
  full_join(c_prp_ipd) %>% 
  select(trial_id, source, arm_n, arm_n_attr, trttype)

# Remove arms with no attrition
d_no_attr <- c_add_ipd %>% filter(arm_n_attr != 0)

# 348 trials
# 262 agg, 86 ipd
d_no_attr %>% count(trial_id, source) %>% count(source)
d_no_attr %>% reframe(across(everything(), ~ sum(is.na(.))))

# Save (remove ipd when trials are out)
write_csv(d_no_attr, "processed_data/df_counts_arms.csv")
