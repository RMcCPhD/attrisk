
# Pt4 of agg bl extraction
# Checking extraction from 2022 update
# No additional data

source("scripts/00_config.R")
source("scripts/00_paths.R")
library(tidyverse)

# Helper: Check missingness
checkNA <- function(df) {
  df %>% 
    reframe(across(everything(), ~ sum(is.na(.)))) %>% 
    pivot_longer(everything()) %>% 
    filter(value > 0)
}

# Import trial set
a_imp_set <- readRDS("processed_data/latest_dataset.rds")
checkNA(a_imp_set)

# Import read_tabletidier dataset
# No additional data
a_imp_pt4 <- read_csv(z_imp_tt) %>% inner_join(a_imp_set %>% distinct(trial_id))
checkNA(a_imp_pt4)
