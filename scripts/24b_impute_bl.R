
# Impute missing variable parameters using amelia
source("scripts/00_config.R")
library(tidyverse)
library(Amelia)

# Import bl chars
a_imp_bl <- readRDS("processed_data/bl_chars_noimp.rds")

# Check missingness
a_imp_bl %>% 
  reframe(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  mutate(prop = value/344) %>% 
  print(n = 25)

# Setup for amelia -------------------------------------------------------------

# Format trttype
b_format <- a_imp_bl %>% mutate(trttype = factor(trttype))

# Remove reference variables - white, former_never
# Will add them as remainder of other imputed variables
b_rm <- b_format %>% 
  select(-c(white, former_never)) %>% 
  mutate(
    across(asian:current, ~ if_else(. == 0, NA, .)),
    across(black:other, ~ if_else(asian > 0.9999, 0, .))
  )

# Log-transformations (so SDs are always positive)
logs <- c(
  "sd_age", "sd_sbp"," sd_bmi", "sd_t2dur", "sd_egfr",
  "sd_hba1c", "sd_hdl", "sd_total_chol"
)

# Set bounds for proportion imputations
prop_vars <- intersect(
  c("female_prop","asian","black","other","current"), 
  names(b_rm)
)

bounds <- do.call(
  rbind, 
  lapply(prop_vars, function(v) c(match(v, names(b_rm)), 0, 1))
)

# Set bounds for mean imputations
mean_limits <- rbind(
  c(match("mean_age", names(b_rm)),        40, 75),
  c(match("mean_sbp", names(b_rm)),        124, 152),
  c(match("mean_bmi", names(b_rm)),        24, 37),
  c(match("mean_t2dur", names(b_rm)),      0.7, 21),
  c(match("mean_egfr", names(b_rm)),       24, 154),
  c(match("mean_hba1c", names(b_rm)),      7, 10),
  c(match("mean_hdl", names(b_rm)),        0.9, 3),
  c(match("mean_total_chol", names(b_rm)), 4, 10)
)

bounds <- rbind(
  bounds, 
  mean_limits[!is.na(mean_limits[,1]), , drop = FALSE]
)

# Run amelia -------------------------------------------------------------------

set.seed(123)

b_impute <- amelia(
  b_rm,
  m = 20,
  idvars = "trial_id",
  noms   = "trttype",
  logs   = intersect(logs, names(b_rm)),
  bounds = bounds,
  empri  = 0.01 * nrow(b_rm),
  incheck = TRUE
)

# Adding back in reference group (if present), compute based on other groups
b_rebuild <- lapply(b_impute$imputations, function(df) {
  b_format %>% 
    select(trial_id, trttype, white, former_never) %>% 
    left_join(df) %>% 
    select(trial_id, trttype, mean_age:other, white, current, former_never) %>% 
    mutate(
      white = if_else(is.na(white), 1 - (asian + black + other), white),
      former_never = if_else(is.na(former_never), 1 - current, former_never)
    )
})

# Take first imputed dataset
# Sanity check on race and smoking status
b_imputed <- b_rebuild$imp1 %>% 
  mutate(
    across(mean_age:former_never, ~ round(., 2)),
    # asian = if_else(asian == 0.01, 1, asian),
    race_total = asian + black + other + white,
    smoke_total = current + former_never
  )
