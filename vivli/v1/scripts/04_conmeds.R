
source("Scripts/00_config.R")
library(tidyverse)
library(haven)

# Import data, dictionary and lookups
a_imp_grps <- readRDS("Processed_data/n95_t2e_fixes.rds") %>% distinct(nct_id, id)
a_imp_dd <- read_rds("../Lookups/Scratch_data/data_dictionary_231.rds")
a_imp_lkp <- read_csv("Created_metadata/lookup_trials.csv")
a_imp_atc <- read_csv("Created_metadata/atc_all_levels.csv") %>% filter(!grepl("^V", atccd4))

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

# Sponsors and trial counts
b_dd_n95 %>% count(sponsor, nct_id) %>% count(sponsor)

# Astrazeneca ------------------------------------------------------------------

# 4 trials in dd, manually extract 5th
b_dd_n95 %>% filter(sponsor == "astrazeneca") %>% count(nct_id)
c_dd <- b_dd_n95 %>% filter(sponsor == "astrazeneca" & grepl("cm", dataset))
c_dd %>% count(nct_id, trial)
c_paths <- c_dd %>% select(nct_id, dataset, filepath) %>% distinct()

# Import datasets
c_lst_imp <- lapply(1:nrow(c_paths), function(i) {
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

# Find common variables among dfs
c_cmn_vars <- Reduce(intersect, lapply(c_lst_imp, names))

# Get common variables
c_cmn_df <- lapply(seq_along(c_lst_imp), function(i) {
  c_lst_imp[[i]] %>% 
    select(c(
      all_of(c_cmn_vars), any_of(c("cmendnfl", "cmenrtpt")), 
      ends_with("atc"), starts_with("atc")
    ))
})

# Join covariate set participants here to simplify extraction
c_cmn_df <- bind_rows(c_cmn_df) %>% 
  rename(
    study = studyid,
    id = usubjid
  ) %>% 
  inner_join(b_grps %>% filter(sponsor == "astrazeneca"))

c_tidy_df <- c_cmn_df %>% 
  select(study, id, cmdecod, contains("atc"), cmendnfl, cmenrtpt) %>% 
  mutate(
    cm = if_else(!is.na(cmdecod), cmdecod, atcdtxt),
    atc = if_else(!is.na(cmspcatc), cmspcatc, atccd),
    status = if_else(!is.na(cmendnfl), cmendnfl, cmenrtpt)
  ) %>% 
  select(study, id, cm, status) %>% 
  filter(status %in% c("0", "ONGOING")) %>% 
  left_join(a_imp_atc %>% select(cm = cmterm, atccd4, atctxt4)) %>% 
  group_by(study, id, cm) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  group_by(study, id, atccd4) %>% 
  reframe(n_atc = n()) %>%
  filter(!is.na(atccd4)) %>% 
  group_by(study, id) %>% 
  reframe(
    unq_atc = paste0(unique(atccd4), collapse = ", "),
    n_atc = length(unique(atccd4))
  )

summary(c_tidy_df$n_atc)

# Manually extract NCT01619059
c_manual_path <- "V:/From 8697/CV181-168/IPD.zip"
c_manual_df <- "IPD/deid_adcm.sas7bdat"

c_tidy_manual <- read_sas(unz(c_manual_path, c_manual_df)) %>% 
  rename_with(., tolower) %>% 
  filter(cmendnfl == 0) %>%
  select(study = studyid, id = usubjid, cm = cmdecod) %>% 
  left_join(a_imp_atc %>% select(cm = cmterm, atccd4, atctxt4)) %>% 
  filter(!is.na(atccd4)) %>% 
  group_by(study, id, cm) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  group_by(study, id, atccd4) %>% 
  reframe(n_atc = n()) %>% 
  group_by(study, id) %>% 
  reframe(
    unq_atc = paste0(unique(atccd4), collapse = ", "),
    n_atc = length(unique(atccd4))
  )

summary(c_tidy_manual$n_atc)

c_tidy <- c_tidy_df %>% full_join(c_tidy_manual)

c_tidy_join <- b_grps %>% 
  filter(sponsor == "astrazeneca") %>% 
  left_join(c_tidy)

# Check
c_tidy_join %>%
  group_by(study, nct_id) %>% 
  reframe(
    n = n(),
    nas = sum(is.na(n_atc)),
    median_cm = median(n_atc, na.rm = TRUE)
  )

# Save
saveRDS(c_tidy_join, "Processed_data/conmeds/az_tidy.rds")

# Remove to save memory
rm(list = ls(pattern = "^c_"))
gc()

# Boehringer -------------------------------------------------------------------

b_dd_n95 %>% filter(sponsor == "boehringer") %>% count(nct_id, trial) %>% print(n = 33)
d_dd <- b_dd_n95 %>% filter(sponsor == "boehringer" & grepl("basco\\.|bascoad\\.", dataset))
d_dd %>% count(nct_id, trial)
d_paths <- d_dd %>% select(nct_id, dataset, filepath) %>% distinct()

# Handle bascoad datasets
d_lst_bascoad <- lapply(1:nrow(d_paths[grepl("bascoad", d_paths$dataset),]), function(i) {
  df <- d_paths[grepl("bascoad", d_paths$dataset),]
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
names(d_lst_bascoad) <- d_paths$dataset[grepl("bascoad", d_paths$dataset)]

d_cmn_vars <- Reduce(intersect, lapply(d_lst_bascoad, names))

d_cmn_df <- lapply(seq_along(d_lst_bascoad), function(i) {
  d_lst_bascoad[[i]] %>% 
    select(c(
      all_of(d_cmn_vars)
    ))
})

d_cmn_df <- bind_rows(d_cmn_df)

d_tidy_df <- d_cmn_df %>% 
  select(study = origin, id = adp_pid, med = baddc) %>% 
  filter(med != "No treatment") %>% 
  separate(
    med, 
    into = c("med1", "med2", "med3", "med4"), 
    sep = " & "
  ) %>% 
  pivot_longer(
    contains("med"),
    names_to = "med_name",
    values_to = "med"
  ) %>% 
  filter(!is.na(med)) %>% 
  select(-med_name) %>% 
  mutate(
    med = if_else(grepl("only", med), gsub(" only", "", med), med),
    atc = case_when(
      med == "A-glucosidase inhib" ~ "A10BF",
      med == "DPP-IV inhibitor" ~ "A10BH",
      med == "GLP-1 agonist" ~ "A10BJ",
      med == "Glinide" ~ "A10BX",
      med == "Glitazone" ~ "A10BG",
      med == "Insulin" ~ "A10A",
      med == "Metformin" ~ "A10BA",
      med == "Other antidiabetics" ~ "A10BX",
      med == "Sulfonylurea" ~ "A10BB"
    )
  ) %>% 
  group_by(study, id, atc) %>% 
  reframe(n_atc = n()) %>% 
  group_by(study, id) %>% 
  reframe(
    unq_atc = paste0(unique(atc), collapse = ", "),
    n_atc = length(unique(atc))
  )

summary(d_tidy_df$n_atc)

# Handle basco datasets
d_lst_basco <- lapply(1:nrow(d_paths[grepl("basco\\.", d_paths$dataset),]), function(i) {
  df <- d_paths[grepl("basco\\.", d_paths$dataset),]
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

d_cmn_vars2 <- Reduce(intersect, lapply(d_lst_basco, names))

d_cmn_df2 <- lapply(seq_along(d_lst_basco), function(i) {
  d_lst_basco[[i]] %>% 
    select(c(
      all_of(d_cmn_vars2), contains("adm")
    ))
})

d_cmn_df2 <- bind_rows(d_cmn_df2)

d_tidy_df2 <- d_cmn_df2  %>%
  select(
    study, id = adp_pid, 
    badm1dc, badm2dc, badm3dc, badm4dc, badm5dc, 
    badm6dc, badm7dc, badm8dc, badm9dc, badm99dc, sadm1dc, 
    sadm2dc, sadm4dc, sadm5dc, sadm6dc, sadm7dc, 
    sadm8dc, sadm99dc
  ) %>% 
  pivot_longer(
    badm1dc:sadm99dc,
    names_to = "med",
    values_to = "med_name"
  ) %>% 
  filter(!(med_name %in% c(
    NA, "", "                   ."
  ))) %>%  
  mutate(
    atc = case_when(
      med_name == "A-glucosidase inhib" | grepl("Alpha", med_name) ~ "A10BF",
      med_name == "DPP-IV inhibitor" ~ "A10BH",
      med_name == "GLP-1 agonist" ~ "A10BJ",
      med_name == "Glinide" ~ "A10BX",
      med_name == "Glitazone" ~ "A10BG",
      med_name == "Insulin" ~ "A10A",
      med_name == "Metformin" ~ "A10BA",
      grepl("Other", med_name) ~ "A10BX",
      grepl("Sulf", med_name) ~ "A10BB"
    )
  ) %>% 
  group_by(study, id, atc) %>% 
  reframe(n_atc = n()) %>% 
  group_by(study, id) %>% 
  reframe(
    unq_atc = paste0(unique(atc), collapse = ", "),
    n_atc = length(unique(atc))
  )

summary(d_tidy_df2$n_atc)

d_tidy <- d_tidy_df %>% full_join(d_tidy_df2)

d_tidy_join <- b_grps %>% 
  filter(sponsor == "boehringer") %>% 
  left_join(d_tidy)

# Check
d_tidy_join %>%
  group_by(study, nct_id) %>% 
  reframe(
    n = n(),
    nas = sum(is.na(n_atc)),
    median_cm = median(n_atc, na.rm = TRUE)
  ) %>% 
  arrange(nct_id) %>% 
  print(n = 33)

# Save
saveRDS(d_tidy_join, "Processed_data/conmeds/boe_tidy.rds")

# Remove to save memory
rm(list = ls(pattern = "^d_"))
gc()

# Eli lilly --------------------------------------------------------------------

b_dd_n95 %>% filter(sponsor == "eli lilly") %>% count(nct_id, trial)

e_dd <- b_dd_n95 %>% 
  filter(
    (sponsor == "eli lilly" & grepl("cmtpyatc|cmtpy|cm\\.", dataset))
    | trial == "13439" & grepl("cmtpy", dataset),
    !grepl("supp", dataset)
  )

e_dd %>% count(nct_id, trial)
e_paths <- e_dd %>% select(nct_id, dataset, filepath) %>% distinct()

# Handle cmptpy datasets
e_lst_cmtpy <- lapply(1:nrow(e_paths[grepl("cmtpy\\.", e_paths$dataset),]), function(i) {
  df <- e_paths[grepl("cmtpy\\.", e_paths$dataset),]
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

e_cmn_vars <- Reduce(intersect, lapply(e_lst_cmtpy, names))

e_cmn_df <- lapply(seq_along(e_lst_cmtpy), function(i) {
  e_lst_cmtpy[[i]] %>% 
    select(c(
      all_of(e_cmn_vars)
    ))
})

e_cmn_df <- bind_rows(e_cmn_df)

e_tidy_df <- e_cmn_df %>% 
  select(study = sdyid, id = usubjid, cmendtc, cmterm) %>% 
  filter(is.na(cmendtc) | cmendtc == "") %>% 
  left_join(a_imp_atc %>% select(cmterm, atccd4, atctxt4)) %>% 
  group_by(study, id, cmterm) %>% 
  slice_head(n = 1) %>% 
  group_by(study, id, atccd4) %>% 
  reframe(n_atc = n()) %>% 
  group_by(study, id) %>% 
  reframe(
    unq_atc = paste0(unique(atccd4), collapse = ", "),
    n_atc = length(unique(atccd4))
  )

summary(e_tidy_df$n_atc)

# Handle cm datasets
e_lst_cm <- lapply(1:nrow(e_paths[!grepl("cmtpy", e_paths$dataset),]), function(i) {
  df <- e_paths[!grepl("cmtpy", e_paths$dataset),]
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

e_cmn_vars3 <- Reduce(intersect, lapply(e_lst_cm, names))

e_cmn_cm <- lapply(seq_along(e_lst_cm), function(i) {
  e_lst_cm[[i]] %>% 
    select(c(
      all_of(e_cmn_vars3)
    ))
})

e_cmn_cm <- bind_rows(e_cmn_cm)

e_tidy_cm <- e_cmn_cm %>% 
  select(study = studyid, id = usubjid, cmdecod, cmendtc) %>% 
  filter(cmendtc == "" | is.na(cmendtc)) %>% 
  left_join(a_imp_atc %>% select(cmdecod = cmterm, atccd4, atctxt4)) %>% 
  group_by(study, id, cmdecod) %>% 
  slice_head(n = 1) %>% 
  group_by(study, id, atccd4) %>% 
  reframe(n_atc = n()) %>% 
  group_by(study, id) %>% 
  reframe(
    unq_atc = paste0(unique(atccd4), collapse = ", "),
    n_atc = length(unique(atccd4))
  )

summary(e_tidy_cm$n_atc)

e_tidy <- e_tidy_df %>% full_join(e_tidy_cm)

e_tidy_join <- b_grps %>% 
  filter(sponsor == "eli lilly") %>% 
  left_join(e_tidy)

# Check
e_tidy_join %>%
  group_by(study, nct_id) %>% 
  reframe(
    n = n(),
    nas = sum(is.na(n_atc)),
    median_cm = median(n_atc, na.rm = TRUE)
  ) %>% 
  arrange(nct_id)

# Save
saveRDS(e_tidy_join, "Processed_data/conmeds/eli_tidy.rds")

# Remove to save memory
rm(list = ls(pattern = "^e_"))
gc()

# GSK --------------------------------------------------------------------------

b_dd_n95 %>% filter(sponsor == "gsk") %>% count(nct_id, trial)

f_dd <- b_dd_n95 %>% 
  filter(
    sponsor == "gsk", 
    grepl("cm", dataset),
    !grepl("y2", dataset)
  )

f_dd %>% count(nct_id, trial)
f_paths <- f_dd %>% select(nct_id, dataset, filepath) %>% distinct()

f_lst_cm <- lapply(1:nrow(f_paths), function(i) {
  filepath <- f_paths$filepath[i]
  dataset <- f_paths$dataset[i]
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

f_cmn_vars <- Reduce(intersect, lapply(f_lst_cm, names))

f_cmn_cm <- lapply(seq_along(f_lst_cm), function(i) {
  f_lst_cm[[i]] %>% 
    select(c(
      all_of(f_cmn_vars)
    ))
})

f_cmn_cm <- bind_rows(f_cmn_cm)

f_tidy_cm <- f_cmn_cm %>% 
  select(study = studyid, id = usubjid, cmdecod, cmendtc) %>% 
  filter(
    cmendtc == "" | is.na(cmendtc),
    cmdecod != "" | !is.na(cmdecod),
    cmdecod != ""
  ) %>% 
  group_by(study, id) %>% 
  separate_rows(cmdecod, sep = ";") %>% 
  ungroup() %>% 
  left_join(a_imp_atc %>% select(cmdecod = cmterm, atccd4, atctxt4)) %>% 
  group_by(study, id, cmdecod) %>% 
  slice_head(n = 1) %>% 
  group_by(study, id, atccd4) %>% 
  reframe(n_atc = n()) %>% 
  group_by(study, id) %>% 
  reframe(
    unq_atc = paste0(unique(atccd4), collapse = ", "),
    n_atc = length(unique(atccd4))
  )

summary(f_tidy_cm$n_atc)

f_tidy_join <- b_grps %>% 
  filter(sponsor == "gsk") %>% 
  left_join(f_tidy_cm)

# Check
f_tidy_join %>%
  group_by(study, nct_id) %>% 
  reframe(
    n = n(),
    nas = sum(is.na(n_atc)),
    median_cm = median(n_atc, na.rm = TRUE)
  ) %>% 
  arrange(nct_id)

# Save
saveRDS(f_tidy_join, "Processed_data/conmeds/gsk_tidy.rds")

# Remove to save memory
rm(list = ls(pattern = "^f_"))
gc()

# Janssen ----------------------------------------------------------------------

b_dd_n95 %>% filter(sponsor == "janssen") %>% count(nct_id, trial)

g_dd <- b_dd_n95 %>% 
  filter(
    sponsor == "janssen",
    grepl("cm|CM", dataset),
    !grepl("supp|SUPP", dataset)
  )

g_dd %>% count(nct_id, trial)
g_paths <- g_dd %>% select(nct_id, dataset, filepath) %>% distinct()

g_lst_imp_cm <- lapply(1:nrow(g_paths), function(i) {
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

g_cmn_vars <- Reduce(intersect, lapply(g_lst_imp_cm, names))
g_cmn_cm <- lapply(seq_along(g_lst_imp_cm), function(i) {
  g_lst_imp_cm[[i]] %>% 
    select(c(any_of(g_cmn_vars)), contains(c("end", "subj")))
})

g_cmn_cm <- bind_rows(g_cmn_cm)

g_tidy_df <- g_cmn_cm %>% 
  mutate(
    usubjid = if_else(!is.na(usubjid), usubjid, dusubjid),
    end = case_when(
      !is.na(cmendtc) & is.na(cmendy) ~ cmendtc,
      is.na(cmendtc) & !is.na(cmendy) ~ cmendy,
      TRUE ~ cmendtc
    )
  ) %>% 
  select(study = studyid, id = usubjid, cmdecod, atc = cmclascd, end) %>% 
  filter(is.na(end) | end == "") %>% 
  left_join(a_imp_atc %>% select(cmdecod = cmterm, atccd4, atctxt4)) %>% 
  group_by(study, id, cmdecod) %>% 
  slice_head(n = 1) %>% 
  group_by(study, id, atccd4) %>% 
  reframe(n_atc = n()) %>% 
  group_by(study, id) %>% 
  reframe(
    unq_atc = paste0(unique(atccd4), collapse = ", "),
    n_atc = length(unique(atccd4))
  ) %>% 
  filter(study != "28431754DIA3009")

summary(g_tidy_df$n_atc)

# Manually get 28431754DIA3009
g_dia3009_path <- "V:/CR016480/cm.xpt"
g_dia3009_cm <- read_xpt(g_dia3009_path) %>% 
  rename_with(., tolower) %>% 
  select(study = studyid, id = usubjid, atc = cmclascd, cmendtc) %>% 
  filter(cmendtc == "") %>% 
  left_join(a_imp_atc %>% select(cmdecod = cmterm, atc = atccd4, atctxt4)) %>% 
  group_by(study, id, cmdecod) %>% 
  slice_head(n = 1) %>% 
  group_by(study, id, atc) %>% 
  reframe(n_atc = n()) %>% 
  group_by(study, id) %>% 
  reframe(
    unq_atc = paste0(unique(atc), collapse = ", "),
    n_atc = length(unique(atc))
  )

summary(g_dia3009_cm$n_atc)

g_tidy_join <- b_grps %>% 
  filter(sponsor == "janssen") %>% 
  left_join(g_tidy_df %>% full_join(g_dia3009_cm))

# Check
g_tidy_join %>%
  group_by(study, nct_id) %>% 
  reframe(
    n = n(),
    nas = sum(is.na(n_atc)),
    median_cm = median(n_atc, na.rm = TRUE)
  ) %>% 
  arrange(nct_id)

# Save
saveRDS(g_tidy_join, "Processed_data/conmeds/janssen_tidy.rds")

# Remove to save memory
rm(list = ls(pattern = "^g_"))
gc()

# Sanofi -----------------------------------------------------------------------

b_dd_n95 %>% filter(sponsor == "sanofi") %>% count(nct_id, trial)

h_dd <- b_dd_n95 %>% 
  filter(
    sponsor == "sanofi",
    grepl("cm|CM", dataset),
    !grepl("supp|SUPP", dataset)
  )

h_dd %>% count(nct_id, trial)
h_paths <- h_dd %>% select(nct_id, dataset, filepath) %>% distinct()

h_lst_imp_cm <- lapply(1:nrow(h_paths), function(i) {
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
  } else {
    read_sas(filepath) %>%
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  }
})

h_cmn_vars <- Reduce(intersect, lapply(h_lst_imp_cm, names))
h_cmn_cm <- lapply(seq_along(h_lst_imp_cm), function(i) {
  h_lst_imp_cm[[i]] %>% 
    select(c(any_of(h_cmn_vars)), contains("atc"))
})

h_cmn_cm <- bind_rows(h_cmn_cm)

h_tidy_df <- h_cmn_cm %>% 
  filter(is.na(cmendy)) %>% 
  select(study = studyid, id = usubjid, cmdecod) %>% 
  left_join(a_imp_atc %>% select(cmdecod = cmterm, atc = atccd4, atctxt4)) %>% 
  group_by(study, id, cmdecod) %>% 
  slice_head(n = 1) %>% 
  group_by(study, id, atc) %>% 
  reframe(n_atc = n()) %>% 
  group_by(study, id) %>% 
  reframe(
    unq_atc = paste0(unique(atc), collapse = ", "),
    n_atc = length(unique(atc))
  )

summary(h_tidy_df$n_atc)

h_tidy_join <- b_grps %>% 
  filter(sponsor == "sanofi") %>% 
  left_join(h_tidy_df)

# Check
h_tidy_join %>%
  group_by(study, nct_id) %>% 
  reframe(
    n = n(),
    nas = sum(is.na(n_atc)),
    median_cm = median(n_atc, na.rm = TRUE)
  ) %>% 
  arrange(nct_id)

# Save
saveRDS(h_tidy_join, "Processed_data/conmeds/sanofi_tidy.rds")

# Remove to save memory
rm(list = ls(pattern = "^h_"))
gc()

# Takeda -----------------------------------------------------------------------

b_dd_n95 %>% filter(sponsor == "takeda") %>% count(nct_id, trial)

i_dd <- b_dd_n95 %>% 
  filter(
    sponsor == "takeda",
    grepl("cm|CM|conmed", dataset),
    !grepl("supp|SUPP", dataset)
  )

i_dd %>% count(nct_id, trial)
# i_paths <- i_dd %>% select(nct_id, dataset, filepath) %>% distinct()
i_paths <- read_csv("Scratch_data/mh_cm/tka_flag_paths_cm.csv")

i_lst_imp <- lapply(1:nrow(i_paths), function(i) {
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
  } else {
    read_sas(filepath) %>%
      rename_with(., tolower) %>% 
      select(where(~!all(is.na(.)))) %>% 
      mutate(across(everything(), as.character))
  }
})
names(i_lst_imp) <- paste0(i_paths$nct_id, "_", i_paths$dataset)

# Tidy cm datasets
i_cmn_vars_cm <- Reduce(intersect, lapply(i_lst_imp[grepl("cm", names(i_lst_imp))], names))
i_cmn_cm <- lapply(seq_along(i_lst_imp[grepl("cm", names(i_lst_imp))]), function(i) {
  df <- i_lst_imp[grepl("cm", names(i_lst_imp))][[i]]
  df %>% 
    select(c(any_of(i_cmn_vars_cm)))
})

i_cmn_cm <- bind_rows(i_cmn_cm)

i_tidy_cm <- i_cmn_cm %>% 
  select(study, id = pt, cmongo, cmprfn) %>% 
  filter(cmongo == "1") %>% 
  select(-cmongo) %>% 
  left_join(a_imp_atc %>% select(cmprfn = cmterm, atc = atccd4, atctxt4)) %>% 
  group_by(study, id, cmprfn) %>% 
  slice_head(n = 1) %>% 
  group_by(study, id, atc) %>% 
  reframe(n_atc = n()) %>% 
  group_by(study, id) %>% 
  reframe(
    unq_atc = paste0(unique(atc), collapse = ", "),
    n_atc = length(unique(atc))
  )

summary(i_tidy_cm$n_atc)

# Tidy d_conmed datasets
i_cmn_vars_conmed <- Reduce(intersect, lapply(i_lst_imp[grepl("d_conmed", names(i_lst_imp))], names))
i_cmn_conmed <- lapply(seq_along(i_lst_imp[grepl("d_conmed", names(i_lst_imp))]), function(i) {
  df <- i_lst_imp[grepl("d_conmed", names(i_lst_imp))][[i]]
  df %>% 
    select(c(
      any_of(i_cmn_vars_conmed)
    ))
})

i_cmn_conmed <- bind_rows(i_cmn_conmed)

i_tidy_conmed <- i_cmn_conmed %>% 
  select(study, id = subjid, cmongo, prefnam) %>% 
  filter(cmongo == "1") %>% 
  select(-cmongo) %>% 
  left_join(a_imp_atc %>% select(prefnam = cmterm, atc = atccd4, atctxt4)) %>% 
  group_by(study, id, prefnam) %>% 
  slice_head(n = 1) %>% 
  group_by(study, id, atc) %>% 
  reframe(n_atc = n()) %>% 
  group_by(study, id) %>% 
  reframe(
    unq_atc = paste0(unique(atc), collapse = ", "),
    n_atc = length(unique(atc))
  )

summary(i_tidy_conmed$n_atc)

# Last one
i_tidy_conmed2 <- i_lst_imp[["NCT01632007_conmed_1_0.sas7bdat"]] %>% 
  select(study, id = subjid, cmongo, prefnam) %>% 
  filter(cmongo == "1") %>% 
  select(-cmongo) %>% 
  left_join(a_imp_atc %>% select(prefnam = cmterm, atc = atccd4, atctxt4)) %>% 
  group_by(study, id, prefnam) %>% 
  slice_head(n = 1) %>% 
  group_by(study, id, atc) %>% 
  reframe(n_atc = n()) %>% 
  group_by(study, id) %>% 
  reframe(
    unq_atc = paste0(unique(atc), collapse = ", "),
    n_atc = length(unique(atc))
  )

summary(i_tidy_conmed2$n_atc)

i_tidy <- i_tidy_cm %>% 
  full_join(i_tidy_conmed) %>% 
  full_join(i_tidy_conmed2)

i_tidy_join <- b_grps %>% 
  filter(sponsor == "takeda") %>% 
  left_join(i_tidy)

# Check
i_tidy_join %>%
  group_by(study, nct_id) %>% 
  reframe(
    n = n(),
    nas = sum(is.na(n_atc)),
    median_cm = median(n_atc, na.rm = TRUE)
  ) %>% 
  arrange(nct_id)

# Save
saveRDS(i_tidy_join, "Processed_data/conmeds/takeda_tidy.rds")

# Remove to save memory
rm(list = ls(pattern = "^i_"))
gc()
