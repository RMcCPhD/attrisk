
source("Scripts/00_config.R")
library(tidyverse)

# Delphi lookup ----------------------------------------------------------------

# Import
a_imp_delphi <- read_csv("Created_metadata/lookup_delphi_multimorb.csv")
a_imp_terms <- read_csv("Created_metadata/meddra_codes_for_mm.csv")

# Adjust names in line with meddra lookup
b_delphi_cond <- a_imp_delphi %>% 
  mutate(
    condition = case_match(
      condition,
      "alcoholism" ~ "drug or alcohol misuse",
      "benign cerebral tumour" ~ "benign cerebral tumours that ca",
      "hearing impairment" ~ "hearing impairment that cannot",
      "musculoskeletal problems due to injury" ~ "long term msk problem from inju",
      "vision impairment" ~ "vision impairment that cannot b",
      .default = condition
    )
  ) %>% 
  distinct()

# Add preferred term variables
b_add_meddra <- b_delphi_cond %>% 
  left_join(a_imp_terms) %>% 
  rename(
    mhterm = `preferred term`,
    pt_code = `pt code`
  ) %>% 
  mutate(mhterm = tolower(mhterm))

# Not totally NA so can remove NAs
b_add_meddra %>% 
  filter(is.na(mhterm)) %>% 
  count(condition)

# Save
write_csv(b_add_meddra %>% filter(!is.na(mhterm)), "Created_metadata/delphi_ai_meddra.csv")

# ATC lookup -------------------------------------------------------------------

# Eli-lilly trial has a series of ATC lookups for general conmeds
# Using that to build an ATC lookup
a_imp_dd <- read_rds("../Lookups/Scratch_data/data_dictionary_231.rds")

# Get cmtpyatc datasets
b_dd_atc <- a_imp_dd %>% 
  filter(grepl("cmtpyatc", dataset)) %>% 
  distinct(filepath, dataset) %>% 
  distinct()

# Handle atc datasets and join
c_imp_lst <- lapply(1:nrow(b_dd_atc), function(i) {
  filepath <- b_dd_atc$filepath[i]
  dataset <- b_dd_atc$dataset[i]
  print(filepath)
  print(dataset)
  if (grepl("\\.zip", filepath)) {
    read_sas(unz(filepath, dataset)) %>%
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  } else {
    read_sas(filepath) %>% 
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  }
})

d_cmn_vars <- Reduce(intersect, lapply(c_imp_lst, names))

d_cmn_atc <- lapply(seq_along(c_imp_lst), function(i) {
  c_imp_lst[[i]] %>% select(all_of(d_cmn_vars))
})

d_cmn_atc <- bind_rows(d_cmn_atc) %>% 
  select(cmterm, whodrgid, contains("atc")) %>% 
  mutate(across(everything(), ~ if_else(. == "", NA, .))) %>% 
  na.omit() %>% 
  distinct()

# # Save as lookup
write_csv(
  d_cmn_atc,
  "Created_metadata/atc_all_levels.csv"
)
