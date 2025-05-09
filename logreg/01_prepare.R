
rm(list = ls())
gc()

library(tidyverse)

# Import
a_imp_df <- read_csv("processed_data/data_both_models.csv")
glimpse(a_imp_df)

# Format, incl treatment type swapped (0 - placebo, 1 - active)
a_imp_fmt <- a_imp_df %>% 
  mutate(
    trttype = if_else(trttype == 1, 0, 1),
    across(male:hba1c, ~ as.numeric(scale(.)) * 0.1),
    across(n_250:reg_em, ~ as.integer(.))
  )

# Split into datasets
b_trialsub <- a_imp_df

# Design-stage -----------------------------------------------------------------

# Remove perfectly collinear variables
b_design <- a_imp_df %>%
  select(-c(male:hba1c, n_abv, arm_100, dur_abv, arm_4, p4, bol, f_abv, pom_abv))

# Vector for formula and covariates (25 covariates)
b_rhs <- paste0(names(b_design[c(9:34)]), collapse = " + ")
b_formula <- as.formula(paste("cbind(arm_n_attr, arm_n_compl) ~", b_rhs))

# Nest by treatment type, fit models
b_nest <- b_design %>% 
  group_by(trttype) %>% 
  nest() %>% 
  mutate(
    mdl = map(data, ~ glm(b_formula, data = .x, family = "binomial"))
  )

# Get outputs
summary(b_nest$mdl[[1]])
summary(b_nest$mdl[[2]])



# Trial submodel ---------------------------------------------------------------