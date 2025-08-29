
# Pt3 of agg bl extraction
# Checking extractions from read_tabletidier

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

# Import first baseline char extraction
# Missingness is two trials not in aact set, manually add at end
a_imp_pt2 <- readRDS("processed_data/agg_bl_second_ext.rds")
checkNA(a_imp_pt2)

# Create a lookup of what each trial is missing
a_na_lkp <- a_imp_pt2 %>%
  pivot_longer(age_mu:total_chol_d) %>%
  filter(is.na(value)) %>%
  distinct() %>%
  mutate(
    var = case_when(
      grepl("age", name) ~ "age",
      grepl("asian|black|other|white", name) ~ "race",
      grepl("bmi", name) ~ "bmi",
      grepl("current|former", name) ~ "smoke",
      grepl("dur_", name) ~ "t2dur",
      grepl("male", name) ~ "sex",
      grepl("sbp", name) ~ "sbp",
      grepl("egfr", name) ~ "egfr",
      grepl("hdl", name) ~ "hdl",
      grepl("chol", name) ~ "total_chol"
    )
  ) %>%
  select(trial_id, var) %>%
  distinct() %>%
  mutate(flag = "1") %>%
  pivot_wider(names_from = var, values_from = flag) %>%
  mutate(
    across(hdl:sex, ~ if_else(is.na(.), "0", .)),
    hba1c = "1"
  )

# write_csv(a_na_lkp, "scratch_data/bl_missing_pt3.csv")

# Import read_tabletidier baseline data
a_imp_tt <- read_csv("../read_tabletidier/data/final_baseline_data.csv")
