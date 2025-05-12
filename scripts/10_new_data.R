
rm(list = ls())
gc()

options(scipen = 999)

library(tidyverse)
library(broom)

# Import
a_imp_df <- read_csv("processed_data/data_both_models.csv")
glimpse(a_imp_df)

# Remove perfectly collinear variables identified in logreg models
# Format variables
# Treatment arm changes to 0 = plc, 1 = ac
b_rem <- a_imp_df %>%
  select(-c(n_abv, arm_100, dur_abv, arm_4, p4, bol, f_abv, pom_abv)) %>% 
  mutate(
    trttype = if_else(trttype == 1, 0, 1),
    across(male:hba1c, ~ as.numeric(scale(.)) * 0.1),
    across(n_250:reg_em, ~ as.integer(.))
  )

glimpse(b_rem)

# Split into design and trial sub-model datasets
names(b_rem)
c_covs <- names(b_rem[c(15:40)]) # 26 covariates

c_design <- list(
  I = nrow(b_rem),
  K = length(c_covs),
  r = b_rem$arm_n_attr,
  N = b_rem$arm_n,
  dt = b_rem$duration,
  arm = as.integer(b_rem$trttype),
  X = as.matrix(b_rem[, c_covs])
)

saveRDS(c_design, "mars/design_model_mk2/stan_dat.rds")

names(b_rem)
d_covs <- names(b_rem[c(10:40)]) # 26 covariates

d_trial <- list(
  I = nrow(b_rem),
  K = length(d_covs),
  r = b_rem$arm_n_attr,
  N = b_rem$arm_n,
  dt = b_rem$duration,
  arm = as.integer(b_rem$trttype),
  X = as.matrix(b_rem[, d_covs])
)

saveRDS(d_trial, "mars/trial_sm_mk2/stan_dat.rds")
