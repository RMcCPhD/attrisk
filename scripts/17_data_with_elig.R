
source("scripts/00_config.R")
library(tidyverse)

# Import datasets
a_imp_df <- readRDS("processed_data/definitive_dataset.rds")
a_imp_elig <- read_csv("processed_data/elig_tidy.csv")

# Address missingness
b_elig_fill <- a_imp_elig %>% 
  mutate(
    age_upper = case_when(
      is.na(age_upper) ~ 150,
      TRUE ~ age_upper
    ),
    bmi_lower = case_when(
      is.na(bmi_lower) ~ ">=20",
      TRUE ~ bmi_lower
    ),
    bmi_upper = case_when(
      is.na(bmi_upper) ~ "<=45",
      TRUE ~ bmi_upper
    ),
    egfr_lower = case_when(
      egfr_lower <= 30 ~ "<=30",
      egfr_lower > 30 ~ ">30",
      is.na(egfr_lower) ~ "None"
    ),
    egfr_upper = case_when(
      egfr_upper <= 60 ~ "<=60",
      egfr_upper > 60 ~ ">60",
      is.na(egfr_upper) ~ "None"
    )
  )

# Categorise variables as binary variables
c_elig_cat <- b_elig_fill %>% 
  mutate(
    # Min age: 0 = 18, 1 = >18
    min_age_cat = case_when(
      age_lower == 18 ~ 0L,
      age_lower > 18 ~ 1L
    ),
    # Max age: 0 = <=80, 1 = >80
    max_age_cat = case_when(
      age_upper <= 80 ~ 0L,
      age_upper > 80 ~ 1L
    ),
    # Min bmi: 0 = Underweight/healthy (<25), 1 = Obese and above (>=25)
    min_bmi_cat = case_when(
      as.numeric(str_extract(bmi_lower, "[0-9]+")) <= 25 ~ 0L,
      TRUE ~ 1L
    ),
    # Max bmi: Obese (<40), 1 = Morbidly obese (>=40)
    max_bmi_cat = case_when(
      as.numeric(str_extract(bmi_upper, "[0-9]+")) <= 40 ~ 0L,
      TRUE ~ 1L
    ),
    # Min egfr: 0 = None, 1 = <=30, 2 = >30
    min_egfr_cat = case_when(
      egfr_lower == "None" ~ 0L,
      egfr_lower == "<=30" ~ 1L,
      TRUE ~ 2L
    ),
    # Max egfr: 0 = None, 1 = <=60, 2 = >60
    max_egfr_cat = case_when(
      egfr_lower == "None" ~ 0L,
      egfr_lower == "<=60" ~ 1L,
      TRUE ~ 2L
    )
  )

# Add row for missing trial (NCT04170998)
c_add_n1 <- c_elig_cat %>% 
  bind_rows(
    tibble(
      trial_id = "NCT04170998",
      hf = 0,
      age_lower = 19,
      age_upper = 150,
      bmi_lower = "18.5",
      bmi_upper = "<40",
      egfr_lower = "None",
      egfr_upper = "None",
      min_age_cat = 1L,
      max_age_cat = 1L,
      min_bmi_cat = 0L,
      max_bmi_cat = 0L,
      min_egfr_cat = 0L,
      max_egfr_cat = 1L
    )
  )

# Join with covariates - 33 covariates (incl participant chars)
# 29 excl. participant chars
d_join <- a_imp_df %>% 
  left_join(
    c_add_n1 %>% 
      select(trial_id, hf, contains("cat")) %>% 
      distinct()
  ) %>% 
  arrange(trial_id) %>% 
  select(trial_id:duration, hf:max_bmi_cat, everything())

# Save
saveRDS(d_join, "processed_data/definitive_dataset_with_elig.rds")
write_csv(d_join, "processed_data/definitive_dataset_with_elig.csv")
