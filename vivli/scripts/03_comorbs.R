
# Index
# a - Imported data
# b - Prepared data dictionary and lookups
# c - Astrazeneca
# d - Boehringer
# e - Eli lilly
# f - GSK
# g - Janssen
# h - Sanofi
# i - Takeda

source("Scripts/00_config.R")
library(tidyverse)
library(haven)

# Import data, dictionary and lookups
a_imp_grps <- readRDS("Processed_data/n95_t2e_fixes.rds") %>% distinct(nct_id, id)
a_imp_dd <- read_rds("../Lookups/Scratch_data/data_dictionary_231.rds")
a_imp_lkp <- read_csv("Created_metadata/lookup_trials.csv")
a_imp_delphi <- read_csv("Created_metadata/delphi_ai_meddra.csv")

# Get data dictionary for 95 trials
b_dd_n95 <- a_imp_grps %>% 
  select(nct_id) %>% 
  distinct() %>% 
  left_join(
    a_imp_lkp %>% 
      select(sponsor, ctgov, trial) %>% 
      rename(nct_id = ctgov) %>% 
      mutate(trial = gsub(",", "", trial))
  ) %>% 
  left_join(a_imp_dd)

# Get groups for filtering to attrition patients
b_grps <- a_imp_grps %>% 
  left_join(
    a_imp_lkp %>% 
      rename(nct_id = ctgov) %>% 
      select(sponsor, nct_id, trial, study) %>% 
      mutate(across(trial:study, ~ gsub(",", "", .)))
  )

# Get delphi flags for adding on
b_delphi_flags <- a_imp_delphi %>% mutate(flag = 1) %>% filter(!is.na(pt_code))

# Sponsors and trial counts
b_dd_n95 %>% count(sponsor, nct_id) %>% count(sponsor)

# Astrazeneca ------------------------------------------------------------------

# 3/5 trials have medical history datasets
c_dd <- b_dd_n95 %>% filter(sponsor == "astrazeneca" & grepl("mh", dataset))
c_dd %>% count(nct_id)
c_paths <- c_dd %>% select(nct_id, dataset, filepath) %>% distinct()

# Import datasets
c_lst_imp <- lapply(seq_along(c_paths), function(i) {
  filepath <- c_paths$filepath[i]
  dataset <- c_paths$dataset[i]
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

# Import variable lookup and get variables to reduce size
c_vars <- read_csv("Created_metadata/vars_mh_cm/az_mh_labels.csv") %>% 
  filter(flag == 1) %>% 
  select(variable, lbl)

c_vars_df <- lapply(seq_along(c_lst_imp), function(i) {
  c_lst_imp[[i]] %>% select(any_of(c_vars$variable), any_of(c("mhenrtpt", "mhptcd")))
})

# Find common variables among dfs
c_cmn_vars <- Reduce(intersect, lapply(c_vars_df, names))

# Get common variables
c_cmn_df <- lapply(seq_along(c_vars_df), function(i) {
  c_vars_df[[i]] %>% select(all_of(c_cmn_vars), any_of(c("mhenrtpt", "mhptcd")))
})
c_cmn_df <- bind_rows(c_cmn_df)

# Tidy into relevant variables
c_tidy_df <- c_cmn_df %>%
  filter(
    mhcat != "DIABETES HISTORY",
    mhoccur %in% c("", "1", "Y"),
    mhenrtpt %in% c("ONGOING", "", NA)
  ) %>% 
  mutate(
    studyid = if_else(studyid == "", "D1689C00014", studyid),
    mhterm = tolower(mhterm),
    mhptcd = as.numeric(mhptcd)
  ) %>% 
  rename(pt_code = mhptcd)

# Handle mhterms and pt codes separately
# If one is there the other tends not to be
c_tidy_mhterm <- c_tidy_df %>% 
  filter(!is.na(mhterm) & mhterm != "") %>% 
  select(-pt_code) %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct()

# Count comorbidities (leave unknowns as NA)
c_tidy_mhpt <- c_tidy_df %>% 
  filter(is.na(mhterm) | mhterm == "") %>% 
  select(-mhterm) %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct()

c_tidy_count <- c_tidy_mhterm %>% 
  full_join(c_tidy_mhpt) %>% 
  filter(!is.na(pivot_name)) %>% 
  select(studyid, usubjid, condition, pivot_name, flag) %>% 
  distinct() %>% 
  group_by(studyid, usubjid, pivot_name) %>%
  mutate(flag = if_else(sum(flag, na.rm = TRUE) >= 1, 1, 0)) %>% 
  distinct() %>% 
  group_by(studyid, usubjid) %>% 
  reframe(comorbs = sum(flag)) %>% 
  rename(
    study = studyid,
    id = usubjid
  )

# Join comorbidity count to attrition patients
c_tidy_join <- b_grps %>% 
  filter(sponsor == "astrazeneca") %>% 
  left_join(c_tidy_count)

# Check
c_tidy_join %>% count(nct_id, comorbs)

# Save
saveRDS(c_tidy_join, "Processed_data/comorbs/az_tidy.rds")

# Remove to save memory
rm(list = ls(pattern = "^c_"))
gc()

# Boehringer -------------------------------------------------------------------

# All have a dataset with medical history or flagged conditions
d_dd <- b_dd_n95 %>% 
  filter(
    sponsor == "boehringer",
    grepl("medh|mh|bcond", dataset) | nct_id == "NCT00915772" & grepl("basco", dataset)
  )

d_dd %>% count(nct_id)

# DDs for each medical history type
d_paths <- d_dd %>% select(nct_id, dataset, filepath) %>% distinct()
# d_paths_flag <- d_paths %>% filter(grepl("bcond|basco", dataset))
d_paths_flag <- read_csv("Scratch_data/mh_cm/boe_flag_paths.csv")
d_paths_medh <- d_paths %>% filter(!grepl("bcond|basco", dataset))

## Boehringer bcond/basco ------------------------------------------------------

# Import datasets
d_lst_imp_flags <- lapply(1:nrow(d_paths_flag), function(i) {
  filepath <- d_paths_flag$filepath[i]
  dataset <- d_paths_flag$dataset[i]
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
names(d_lst_imp_flags) <- d_paths_flag$dataset

# Get relevant variables
d_vars_df <- lapply(seq_along(d_lst_imp_flags), function(i) {
  d_lst_imp_flags[[i]] %>% 
    select(contains(c(
      "study", "adp_pid", 
      "cad", "cvd", "diabft", "diabne", "diabnr", 
      "diabre", "hypert", "paod", "mpt"
    )))
})
names(d_vars_df) <- d_paths_flag$dataset

# Join datasets (19 and 1 trials respectively)
d_bcond <- bind_rows(d_vars_df[grepl("bcond", names(d_vars_df))])
d_basco <- bind_rows(d_vars_df[!grepl("bcond", names(d_vars_df))])

# Tidy bcond, join delphi lookup and count comorbidities
d_tidy_bcond <- d_bcond %>% 
  mutate(
    mhterm = if_else(!is.na(adp_cdlltcd_mpt), adp_cdlltcd_mpt, adp_mpt),
    mhterm = tolower(mhterm),
    pt_code = if_else(!is.na(adp_cdlltcd_mptcd), adp_cdlltcd_mptcd, adp_mptcd),
    pt_code = as.numeric(pt_code)
  ) %>% 
  select(study, id = adp_pid, mhterm, pt_code) %>% 
  filter(mhterm != "" & pt_code != "") %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct() %>% 
  filter(!is.na(flag)) %>% 
  group_by(study, id, pivot_name) %>%
  mutate(flag = if_else(sum(flag, na.rm = TRUE) >= 1, 1, 0)) %>% 
  distinct() %>% 
  group_by(study, id) %>% 
  reframe(comorbs = sum(flag))

# Tidy basco in same manner
d_tidy_basco <- d_basco %>% 
  select(!ends_with("dc"), -c(prestudy, contains("diab"))) %>%
  pivot_longer(
    cad:paod,
    names_to = "condition",
    values_to = "flag"
  ) %>% 
  mutate(
    condition = case_match(
      condition,
      "cad" ~ "coronary artery disease",
      "cvd" ~ "stroke",
      "hypert" ~ "hypertension",
      "paod" ~ "peripheral artery disease"
    )
  ) %>% 
  filter(flag == "1") %>% 
  select(-flag) %>% 
  distinct() %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct() %>% 
  filter(!is.na(flag)) %>% 
  rename(id = adp_pid) %>% 
  select(study, id, pivot_name, flag) %>%
  group_by(study, id, pivot_name) %>% 
  mutate(flag = if_else(sum(flag, na.rm = TRUE) >= 1, 1, 0)) %>% 
  distinct() %>% 
  group_by(study, id) %>% 
  reframe(comorbs = sum(flag))

# Join together
d_tidy1 <- d_tidy_bcond %>% full_join(d_tidy_basco)

## Boehringer medh -------------------------------------------------------------

d_lst_imp_medh <- lapply(1:nrow(d_paths_medh), function(i) {
  filepath <- d_paths_medh$filepath[i]
  dataset <- d_paths_medh$dataset[i]
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
names(d_lst_imp_medh) <- d_paths_medh$dataset

d_cmn_vars <- Reduce(intersect, lapply(d_lst_imp_medh, names))
d_cmn_df <- lapply(seq_along(d_lst_imp_medh), function(i) {
  d_lst_imp_medh[[i]] %>% select(all_of(d_cmn_vars))
})

d_cmn_df <- bind_rows(d_cmn_df)

d_tidy_df <- d_cmn_df %>% 
  select(study, id = adp_pid, artery:hypert) %>% 
  pivot_longer(
    artery:hypert,
    names_to = "condition",
    values_to = "flag"
  ) %>%
  mutate(
    condition = case_match(
      condition,
      "artery" ~ "coronary artery disease",
      "cvdis" ~ "stroke",
      "hypert" ~ "hypertension",
      "paod" ~ "peripheral artery disease"
    )
  ) %>% 
  filter(flag == "1") %>% 
  select(-flag) %>% 
  distinct() %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct() %>% 
  filter(!is.na(flag)) %>% 
  select(study, id, pivot_name, flag) %>%
  group_by(study, id, pivot_name) %>% 
  mutate(flag = if_else(sum(flag, na.rm = TRUE) >= 1, 1, 0)) %>% 
  distinct() %>% 
  group_by(study, id) %>% 
  reframe(comorbs = sum(flag))

## Boehringer trials with flagged and medDRA data ------------------------------

# Identify studies (16)
sum(unique(d_cmn_df$study) %in% unique(d_bcond$study))
d_ids <- d_cmn_df %>% distinct(study) %>% inner_join(d_bcond %>% distinct(study))

# Repeat extraction just for the 16 trials
d_bcond_n16 <- d_bcond %>% 
  inner_join(d_ids) %>% 
  mutate(
    mhterm = if_else(!is.na(adp_cdlltcd_mpt), adp_cdlltcd_mpt, adp_mpt),
    mhterm = tolower(mhterm),
    pt_code = if_else(!is.na(adp_cdlltcd_mptcd), adp_cdlltcd_mptcd, adp_mptcd),
    pt_code = as.numeric(pt_code)
  ) %>% 
  select(study, id = adp_pid, mhterm, pt_code) %>% 
  filter(mhterm != "" & pt_code != "") %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct() %>% 
  filter(!is.na(flag)) %>% 
  group_by(study, id, pivot_name) %>%
  mutate(flag = if_else(sum(flag, na.rm = TRUE) >= 1, 1, 0)) %>% 
  distinct() %>% 
  ungroup() %>% 
  select(study, id, pivot_name, flag)

d_medh_n16 <- d_cmn_df %>% 
  inner_join(d_ids) %>% 
  select(study, id = adp_pid, artery:hypert) %>% 
  pivot_longer(
    artery:hypert,
    names_to = "condition",
    values_to = "flag"
  ) %>%
  mutate(
    condition = case_match(
      condition,
      "artery" ~ "coronary artery disease",
      "cvdis" ~ "stroke",
      "hypert" ~ "hypertension",
      "paod" ~ "peripheral artery disease"
    )
  ) %>% 
  filter(flag == "1") %>% 
  select(-flag) %>% 
  distinct() %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct() %>% 
  filter(!is.na(flag)) %>% 
  select(study, id, pivot_name, flag) %>%
  group_by(study, id, pivot_name) %>% 
  mutate(flag = if_else(sum(flag, na.rm = TRUE) >= 1, 1, 0)) %>% 
  distinct() %>% 
  ungroup()

# Join together and take distinct entries per participant
d_n16_join <- d_bcond_n16 %>% 
  full_join(d_medh_n16) %>% 
  distinct() %>% 
  group_by(study, id) %>% 
  reframe(comorbs = sum(flag)) %>% 
  ungroup()

## Boehringer join -------------------------------------------------------------

# Join all together
d_tidy2 <- d_tidy1 %>% 
  full_join(d_tidy_df %>% filter(!(study %in% d_tidy1$study)))

# Join tidy bcond with groups
d_tidy_join <- b_grps %>% 
  filter(sponsor == "boehringer" & study %in% d_tidy2$study) %>% 
  left_join(d_tidy2) %>% 
  arrange(nct_id, id)

# Check before fix for 16 trials
d_tidy_join %>% 
  count(nct_id, comorbs) %>% 
  ggplot(aes(x = as.integer(comorbs), y = n)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_x_continuous(n.breaks = 10)

# Remove 16 trials from original tidied version
# Add new data
d_tidy_new <- d_tidy_join %>% 
  anti_join(d_ids) %>% 
  full_join(
    d_n16_join %>% 
      left_join(b_grps) %>% 
      group_by(study) %>% 
      fill(nct_id:trial, .direction = "downup") %>% 
      ungroup()
  ) %>% 
  arrange(nct_id, id)

# Check
d_tidy_new %>% 
  count(nct_id, comorbs) %>% 
  ggplot(aes(x = as.integer(comorbs), y = n)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_x_continuous(n.breaks = 10)

# Save
saveRDS(d_tidy_new, "Processed_data/comorbs/boe_tidy_new.rds")

# Remove to save memory
rm(list = ls(pattern = "^d_"))
gc()

# Eli lilly --------------------------------------------------------------------

e_dd <- b_dd_n95 %>% 
  filter(
    sponsor == "eli lilly",
    grepl("history|mh\\.", dataset)
  )

e_dd %>% count(nct_id)

e_paths <- e_dd %>% select(nct_id, dataset, filepath) %>% distinct()
unique(e_paths$dataset)

# History datasets
e_lst_imp_hist <- lapply(1:nrow(e_paths[grepl("history", e_paths$dataset),]), function(i) {
  df <- e_paths[grepl("history", e_paths$dataset),]
  filepath <- df$filepath[i]
  dataset <- df$dataset[i]
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

e_cmn_vars <- Reduce(intersect, lapply(e_lst_imp_hist, names))
e_cmn_hist <- lapply(seq_along(e_lst_imp_hist), function(i) {
  e_lst_imp_hist[[i]] %>% 
    select(any_of(e_cmn_vars))
})

e_cmn_hist <- bind_rows(e_cmn_hist)

e_tidy_hist <- e_cmn_hist %>% 
  select(sdyid, usubjid, pterm, hxendt, hxendtc) %>% 
  filter(is.na(hxendt) | hxendtc == "") %>% 
  select(-hxendt, -hxendtc) %>% 
  mutate(mhterm = tolower(pterm)) %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct() %>% 
  filter(!is.na(flag)) %>% 
  select(study = sdyid, usubjid, pivot_name, flag) %>%
  group_by(study, usubjid, pivot_name) %>% 
  mutate(flag = if_else(sum(flag, na.rm = TRUE) >= 1, 1, 0)) %>% 
  distinct() %>% 
  group_by(study, usubjid) %>% 
  reframe(comorbs = sum(flag))

# Mh datasets
e_lst_imp_mh <- lapply(1:nrow(e_paths[!grepl("history|supp", e_paths$dataset),]), function(i) {
  df <- e_paths[!grepl("history|supp", e_paths$dataset),]
  filepath <- df$filepath[i]
  dataset <- df$dataset[i]
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

e_cmn_vars <- Reduce(intersect, lapply(e_lst_imp_mh, names))
e_cmn_mh <- lapply(seq_along(e_lst_imp_mh), function(i) {
  e_lst_imp_mh[[i]] %>% 
    select(any_of(e_cmn_vars))
})

e_cmn_mh <- bind_rows(e_cmn_mh)

e_tidy_mh <- e_cmn_mh %>% 
  select(studyid, usubjid, mhterm, mhdecod, mhendtc) %>% 
  filter(is.na(mhendtc) | mhendtc == "") %>% 
  select(-mhendtc) %>% 
  mutate(
    mhterm = if_else(mhterm == "" & mhdecod != "", mhdecod, mhterm),
    mhterm = tolower(mhterm)
  ) %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct() %>% 
  filter(!is.na(flag)) %>% 
  select(studyid, usubjid, pivot_name, flag) %>%
  group_by(studyid, usubjid, pivot_name) %>% 
  mutate(flag = if_else(sum(flag, na.rm = TRUE) >= 1, 1, 0)) %>% 
  distinct() %>% 
  group_by(studyid, usubjid) %>% 
  reframe(comorbs = sum(flag))

# Join all together
e_tidy_df <- e_tidy_hist %>% 
  full_join(e_tidy_mh %>% rename(study = studyid)) %>% 
  rename(id = usubjid)

# Join tidy bcond with groups
e_tidy_join <- b_grps %>% 
  filter(sponsor == "eli lilly") %>% 
  left_join(e_tidy_df) %>% 
  arrange(nct_id, id)

# Check
e_tidy_join %>% count(nct_id, comorbs)

# Save
saveRDS(e_tidy_join, "Processed_data/comorbs/eli_tidy.rds")

# Remove to save memory
rm(list = ls(pattern = "^e_"))
gc()

# GSK --------------------------------------------------------------------------

f_dd <- b_dd_n95 %>% 
  filter(
    sponsor == "gsk",
    grepl("mh\\.", dataset),
    !grepl("y2|famh", dataset)
  )

f_dd %>% count(nct_id)

f_paths <- f_dd %>% select(nct_id, dataset, filepath) %>% distinct()
unique(f_paths$dataset)

f_lst_imp_mh <- lapply(1:nrow(f_paths), function(i) {
  df <- f_paths
  filepath <- df$filepath[i]
  dataset <- df$dataset[i]
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
names(f_lst_imp_mh) <- f_paths$dataset

f_cmn_vars <- Reduce(intersect, lapply(f_lst_imp_mh, names))
f_cmn_mh <- lapply(seq_along(f_lst_imp_mh), function(i) {
  f_lst_imp_mh[[i]] %>% 
    select(any_of(f_cmn_vars))
})

f_cmn_mh <- bind_rows(f_cmn_mh)

f_tidy_mh <- f_cmn_mh %>% 
  select(studyid, usubjid, mhterm, mhoccur, mhendtc) %>% 
  filter(mhoccur %in% c("Y", "")) %>% 
  filter(is.na(mhendtc) | mhendtc == "") %>% 
  select(-mhoccur, -mhendtc) %>% 
  mutate(mhterm = tolower(mhterm)) %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct() %>% 
  filter(!is.na(flag)) %>% 
  select(studyid, usubjid, pivot_name, flag) %>%
  group_by(studyid, usubjid, pivot_name) %>% 
  mutate(flag = if_else(sum(flag, na.rm = TRUE) >= 1, 1, 0)) %>% 
  distinct() %>% 
  group_by(studyid, usubjid) %>% 
  reframe(comorbs = sum(flag))

f_tidy_join <- b_grps %>% 
  filter(sponsor == "gsk") %>% 
  left_join(f_tidy_mh %>% rename(study = studyid, id = usubjid)) %>% 
  arrange(nct_id, id)

# Check
f_tidy_join %>% count(nct_id, comorbs)

# Save
saveRDS(f_tidy_join, "Processed_data/comorbs/gsk_tidy.rds")

# Remove to save memory
rm(list = ls(pattern = "^f_"))
gc()

# Janssen ----------------------------------------------------------------------

g_dd <- b_dd_n95 %>% 
  filter(
    sponsor == "janssen",
    grepl("mh\\.|MH\\.", dataset),
    !grepl("supp|SUPP", dataset)
  )

g_dd %>% count(nct_id)

g_paths <- g_dd %>% select(nct_id, dataset, filepath) %>% distinct()
unique(g_paths$dataset)

g_lst_imp_mh <- lapply(1:nrow(g_paths), function(i) {
  df <- g_paths
  filepath <- df$filepath[i]
  dataset <- df$dataset[i]
  print(filepath)
  print(dataset)
  if (grepl("\\.zip", filepath)  & grepl("sas", dataset)) {
    read_sas(unz(filepath, dataset)) %>%
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  } else if (grepl("\\.zip", filepath) & grepl("xpt", dataset)) {
    read_xpt(filepath) %>% 
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  } else if (!grepl("\\.zip", filepath) & grepl("sas", dataset)) {
    read_sas(filepath) %>%
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  } else {
    read_xpt(filepath) %>%
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  }
})

g_cmn_vars <- Reduce(intersect, lapply(g_lst_imp_mh, names))
g_cmn_mh <- lapply(seq_along(g_lst_imp_mh), function(i) {
  g_lst_imp_mh[[i]] %>% 
    select(c(any_of(g_cmn_vars), contains(c("subj", "mhterm"))))
})

g_cmn_mh <- bind_rows(g_cmn_mh)

g_tidy_mh <- g_cmn_mh %>% 
  mutate(
    mhterm = if_else(mhterm %in% c(NA, "") & mhdecod != "", mhdecod, mhterm),
    mhterm = tolower(mhterm),
    usubjid = if_else(is.na(usubjid) & !is.na(dusubjid), dusubjid, usubjid)
  ) %>% 
  select(studyid, usubjid, mhterm, mhoccur) %>% 
  filter(mhoccur %in% c("Y", "")) %>% 
  select(-mhoccur) %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct() %>% 
  filter(!is.na(flag)) %>% 
  select(studyid, usubjid, pivot_name, flag) %>%
  group_by(studyid, usubjid, pivot_name) %>% 
  mutate(flag = if_else(sum(flag, na.rm = TRUE) >= 1, 1, 0)) %>% 
  distinct() %>% 
  group_by(studyid, usubjid) %>% 
  reframe(comorbs = sum(flag)) %>% 
  filter(studyid != "28431754DIA3008")

# Extract 28431754DIA3009 separately, missing from data dictionary
g_dia3009_path <- "V:/CR016480/mh.xpt"

g_dia3009_mh <- read_xpt(g_dia3009_path) %>% 
  rename_with(., tolower) %>% 
  select(studyid, usubjid, mhterm, mhdecod, mhoccur) %>% 
  mutate(
    mhterm = if_else(mhterm == "" & mhdecod != "", mhdecod, mhterm),
    mhterm = tolower(mhterm)
  ) %>% 
  select(studyid, usubjid, mhterm, mhoccur) %>% 
  filter(mhoccur %in% c("Y", "")) %>% 
  select(-mhoccur) %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct() %>% 
  filter(!is.na(flag)) %>% 
  select(studyid, usubjid, pivot_name, flag) %>%
  group_by(studyid, usubjid, pivot_name) %>% 
  mutate(flag = if_else(sum(flag, na.rm = TRUE) >= 1, 1, 0)) %>% 
  distinct() %>% 
  group_by(studyid, usubjid) %>% 
  reframe(comorbs = sum(flag))

# Join together
g_tidy_mh <- g_tidy_mh %>% full_join(g_dia3009_mh)

# Join tidy bcond with groups
g_tidy_join <- b_grps %>% 
  filter(sponsor == "janssen") %>% 
  left_join(g_tidy_mh %>% rename(study = studyid, id = usubjid)) %>% 
  arrange(nct_id, id)

# Check
g_tidy_join %>% count(nct_id, study, comorbs)

# Save
saveRDS(g_tidy_join, "Processed_data/comorbs/janssen_tidy.rds")

# Remove to save memory
rm(list = ls(pattern = "^g_"))
gc()

# Sanofi -----------------------------------------------------------------------

h_dd <- b_dd_n95 %>% 
  filter(
    sponsor == "sanofi",
    grepl("mh\\.", dataset),
    !grepl("supp", dataset)
  )

h_dd %>% count(nct_id)

h_paths <- h_dd %>% select(nct_id, dataset, filepath) %>% distinct()
unique(h_paths$dataset)

h_lst_imp_mh <- lapply(1:nrow(h_paths), function(i) {
  df <- h_paths
  filepath <- df$filepath[i]
  dataset <- df$dataset[i]
  print(filepath)
  print(dataset)
  if (grepl("\\.zip", filepath)  & grepl("sas", dataset)) {
    read_sas(unz(filepath, dataset)) %>%
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  } else if (grepl("\\.zip", filepath) & grepl("xpt", dataset)) {
    read_xpt(filepath) %>% 
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  } else if (!grepl("\\.zip", filepath) & grepl("sas", dataset)) {
    read_sas(filepath) %>%
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  } else {
    read_xpt(filepath) %>%
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  }
})

h_cmn_vars <- Reduce(intersect, lapply(h_lst_imp_mh, names))
h_cmn_mh <- lapply(seq_along(h_lst_imp_mh), function(i) {
  h_lst_imp_mh[[i]] %>% 
    select(c(any_of(h_cmn_vars), contains("mhterm")))
})

h_cmn_mh <- bind_rows(h_cmn_mh)

h_tidy_mh <- h_cmn_mh %>% 
  mutate(
    mhterm = if_else(mhterm %in% c(NA, "") & mhdecod != "", mhdecod, mhterm),
    mhterm = tolower(mhterm)
  ) %>% 
  select(studyid, usubjid, mhterm, mhoccur) %>% 
  filter(mhoccur %in% c("Y", "")) %>% 
  select(-mhoccur) %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct() %>% 
  filter(!is.na(flag)) %>% 
  select(studyid, usubjid, pivot_name, flag) %>%
  group_by(studyid, usubjid, pivot_name) %>% 
  mutate(flag = if_else(sum(flag, na.rm = TRUE) >= 1, 1, 0)) %>% 
  distinct() %>% 
  group_by(studyid, usubjid) %>% 
  reframe(comorbs = sum(flag))

# Join tidy bcond with groups
h_tidy_join <- b_grps %>% 
  filter(sponsor == "sanofi") %>% 
  left_join(h_tidy_mh %>% rename(study = studyid, id = usubjid)) %>% 
  arrange(nct_id, id)

# Check
h_tidy_join %>% count(nct_id)
h_tidy_join %>% count(nct_id, study, comorbs)

# Save
saveRDS(h_tidy_join, "Processed_data/comorbs/sanofi_tidy.rds")

# Remove to save memory
rm(list = ls(pattern = "^h_"))
gc()

# Takeda -----------------------------------------------------------------------

i_dd <- b_dd_n95 %>% 
  filter(
    sponsor == "takeda",
    grepl("mh\\.|medhis", dataset),
    !grepl("supp|a_mh", dataset)
  )

i_dd %>% count(nct_id)

i_paths <- i_dd %>% select(nct_id, dataset, filepath) %>% distinct()
# i_paths <- read_csv("Scratch_data/tka_flag_paths.csv")
unique(i_paths$dataset)

i_lst_imp_mh <- lapply(1:nrow(i_paths), function(i) {
  df <- i_paths
  filepath <- df$filepath[i]
  dataset <- df$dataset[i]
  print(filepath)
  print(dataset)
  if (grepl("\\.zip", filepath)  & grepl("sas", dataset)) {
    read_sas(unz(filepath, dataset)) %>%
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  } else if (grepl("\\.zip", filepath) & grepl("xpt", dataset)) {
    read_xpt(filepath) %>% 
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  } else if (!grepl("\\.zip", filepath) & grepl("sas", dataset)) {
    read_sas(filepath) %>%
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  } else {
    read_xpt(filepath) %>%
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  }
})

i_cmn_vars <- c(
  "study", "pt", "mhongo", "mhprfn", "subjid", "pt_name", "conedty", "conedt",
  "mhendty"
)

i_cmn_mh <- lapply(seq_along(i_lst_imp_mh), function(i) {
  i_lst_imp_mh[[i]] %>% 
    select(any_of(i_cmn_vars))
})

i_cmn_mh <- bind_rows(i_cmn_mh) %>% 
  mutate(
    usubjid = if_else(!is.na(pt) & is.na(subjid), pt, subjid),
    end = case_when(
      !is.na(mhendty) ~ mhendty,
      is.na(mhendty) & !is.na(conedty) ~ conedty,
      is.na(mhendty) & is.na(conedty) & !is.na(conedt) ~ conedt,
      TRUE ~ "CCT_002_ongoing"
    ),
    mhterm = case_when(
      !is.na(mhprfn) & is.na(pt_name) ~ mhprfn,
      TRUE ~ pt_name
    )
  ) %>% 
  select(studyid = study, usubjid, mhterm, end)

i_tidy_mh <- i_cmn_mh %>% 
  filter(end == "" | end == "CCT_002_ongoing") %>% 
  mutate(mhterm = tolower(mhterm)) %>% 
  select(studyid, usubjid, mhterm) %>% 
  left_join(b_delphi_flags %>% filter(condition != "diabetes")) %>% 
  distinct() %>% 
  filter(!is.na(flag)) %>% 
  select(studyid, usubjid, pivot_name, flag) %>%
  group_by(studyid, usubjid, pivot_name) %>% 
  mutate(flag = if_else(sum(flag, na.rm = TRUE) >= 1, 1, 0)) %>% 
  distinct() %>% 
  group_by(studyid, usubjid) %>% 
  reframe(comorbs = sum(flag))

# Join tidy bcond with groups
i_tidy_join <- b_grps %>% 
  filter(sponsor == "takeda") %>% 
  left_join(i_tidy_mh %>% rename(study = studyid, id = usubjid)) %>% 
  arrange(nct_id, id)

# Check
i_tidy_join %>% count(nct_id, study, comorbs)

# Save
saveRDS(i_tidy_join, "Processed_data/comorbs/takeda_tidy.rds")

# Remove to save memory
rm(list = ls(pattern = "^h_"))
gc()
