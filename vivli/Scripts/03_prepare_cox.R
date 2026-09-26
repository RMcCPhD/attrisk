
source("Scripts/00_config.R")
library(tidyverse)
library(survival)
library(broom)

# Import dataset
a_imp_df <- readRDS("Processed_data/02_t2e_checked.rds")

# Add trttype
b_trttype <- a_imp_df %>% 
  mutate(
    trttype = case_when(
      trtcls == "placebo" ~ "plc",
      TRUE ~ "ac"
    )
  ) %>% 
  select(
    sponsor, nct_id, id, arm, trtcls5, trtcls, trttype, 
    everything()
  )

# Check number of events per arm and trial
c_events <- b_trttype %>% 
  group_by(nct_id, trttype) %>% 
  reframe(n_events = sum(event)) %>% 
  mutate(lo = if_else(n_events < 10, 1, 0))

# 15 trials where at least one arm has <10 events
c_events %>% filter(lo == 1) %>% count(nct_id)

# Remove trials with any arm <10 events - 76 trials; 19,023 participants
c_rm_lt10 <- b_trttype %>% 
  anti_join(
    c_events %>%
      filter(lo == 1) %>% 
      distinct(nct_id)
  )

# Add persontime and trial mean persontime (both in years)
d_pt <- c_rm_lt10 %>% 
  group_by(nct_id) %>% 
  mutate(
    pyr = weeks/52,
    mft = sum(pyr)/n()
  ) %>% 
  ungroup() %>% 
  select(
    sponsor, nct_id, id, arm, trtcls5, trtcls, trttype,
    days, weeks, event, category, pyr, mft,
    everything()
  )

d_pt %>% 
  reframe(across(everything(), ~sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  print(n = Inf)

# Save
saveRDS(d_pt, "Processed_data/03_mt10_events_n76.rds")
