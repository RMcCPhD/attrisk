
source("Scripts/00_config.R")
library(tidyverse)
library(haven)

# Import dataset
a_imp_n95 <- readRDS("Processed_data/n95_t2e_fixes.rds")

# Comorbidities ----------------------------------------------------------------

# Import materials
a_imp_files <- list.files("Processed_data/comorbs", full.names = TRUE)
a_imp_paths <- list.files("Processed_data/comorbs")
a_imp_nms <- gsub("_tidy\\.rds", "", a_imp_paths)
a_imp_comorbs <- lapply(a_imp_files, readRDS)
names(a_imp_comorbs) <- a_imp_nms

# Bind data
b_comorbs <- bind_rows(a_imp_comorbs) %>% arrange(sponsor, nct_id, id)

# Add 0 where participants did not have any comorbidities in the set
b_comorbs_zeros <- b_comorbs %>% 
  group_by(nct_id) %>% 
  mutate(
    comorbs = if_else(is.na(comorbs) & sum(is.na(comorbs)) != n(), 0, comorbs)
  ) %>% 
  ungroup()

# Check how many trials don't have comorbidity counts
# 3 don't have medical history
# 3 pt was completely blank
# 3 had no relevant ongoing medical conditions
b_comorbs_zeros %>% 
  group_by(sponsor, nct_id) %>% 
  reframe(nas = sum(is.na(comorbs)) == n()) %>% 
  filter(nas == TRUE)

# Check
summary(b_comorbs_zeros$comorbs)
b_comorbs_zeros %>% ggplot(aes(x = comorbs)) + geom_density()
b_comorbs_zeros %>% ggplot(aes(x = comorbs)) + geom_density() + facet_wrap(~sponsor)
b_comorbs_zeros %>% ggplot(aes(x = comorbs)) + geom_density() + facet_wrap(~nct_id)

# Join with covariates for 92 trials
c_join_comorbs <- a_imp_n95 %>% 
  left_join(b_comorbs_zeros %>% select(nct_id, id, comorbs))

# Recheck in case of id mismatches
c_join_comorbs %>% 
  group_by(sponsor, nct_id) %>% 
  reframe(nas = sum(is.na(comorbs)) == n()) %>% 
  filter(nas == TRUE)

# Save
saveRDS(c_join_comorbs, "Processed_data/n95_comorbs.rds")

# Con meds ---------------------------------------------------------------------

# Import dataset with comorbidities added
a_imp_n95 <- readRDS("Processed_data/n95_comorbs.rds")

# Import materials
a_imp_files <- list.files("Processed_data/conmeds", full.names = TRUE)
a_imp_paths <- list.files("Processed_data/conmeds")
a_imp_nms <- gsub("_tidy\\.rds", "", a_imp_paths)
a_imp_conmeds <- lapply(a_imp_files, readRDS)
names(a_imp_conmeds) <- a_imp_nms

b_conmeds <- bind_rows(a_imp_conmeds) %>% arrange(sponsor, nct_id, id)

b_conmeds_zeros <- b_conmeds %>% 
  group_by(nct_id) %>% 
  mutate(
    n_atc = if_else(is.na(n_atc) & sum(is.na(n_atc)) != n(), 0, n_atc)
  ) %>% 
  ungroup()

# 3 don't have con meds - only treatment naive participants
b_conmeds_zeros %>% 
  group_by(sponsor, nct_id) %>% 
  reframe(nas = sum(is.na(n_atc)) == n()) %>% 
  filter(nas == TRUE)

summary(b_conmeds_zeros$n_atc)
b_conmeds_zeros %>% ggplot(aes(x = n_atc)) + geom_density()
b_conmeds_zeros %>% ggplot(aes(x = n_atc)) + geom_density() + facet_wrap(~sponsor)
b_conmeds_zeros %>% ggplot(aes(x = n_atc)) + geom_density() + facet_wrap(~nct_id)

c_join_conmeds <- a_imp_n95 %>% 
  left_join(b_conmeds_zeros %>% select(nct_id, id, unq_atc, n_atc))

c_join_conmeds %>% 
  group_by(sponsor, nct_id) %>% 
  reframe(nas = sum(is.na(n_atc)) == n()) %>% 
  filter(nas == TRUE)

# Save
saveRDS(c_join_conmeds, "Processed_data/n95_conmeds.rds")
