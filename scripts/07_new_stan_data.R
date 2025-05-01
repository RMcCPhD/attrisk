
# Preparing stan datasets for design model and trial sub-model
# Rerunning design model with changes to arm counts 
# i.e. Combining a few arms across trials found when adding covariates
# Dataset for trial sub-model

source("scripts/00_functions_packages_etc.R")

# Scale participant characteristics
# Mean of 0, SD of 3
a_imp_df <- read_csv("processed_data/data_both_models.csv") %>% 
  arrange(trial_id, trttype) %>% 
  mutate(
    across(n_250:reg_em, as.integer),
    male_prop = as.numeric(scale(male_prop)),
    age = as.numeric(scale(age)),
    bmi = as.numeric(scale(bmi)),
    t2dur = as.numeric(scale(t2dur)),
    hba1c = as.numeric(scale(hba1c))
  )

a_imp_df %>% count(trial_id)
glimpse(a_imp_df)

# Design model dataset
names(a_imp_df)
b_covs <- names(a_imp_df[c(15:48)])

b_design_stan <- list(
  I = nrow(a_imp_df),
  K = length(b_covs),
  r = a_imp_df$arm_n_attr,
  N = a_imp_df$arm_n,
  dt = a_imp_df$duration,
  arm = as.integer(a_imp_df$trttype),
  X = as.matrix(a_imp_df[, b_covs])
)

glimpse(b_design_stan)

# Save
write_rds(b_design_stan, "processed_data/design_stan_data.rds")

# Trial sub-model dataset
names(a_imp_df)
c_covs <- names(a_imp_df[c(10:48)])

c_trial_stan <- list(
  I = nrow(a_imp_df),
  K = length(c_covs),
  r = a_imp_df$arm_n_attr,
  N = a_imp_df$arm_n,
  dt = a_imp_df$duration,
  arm = as.integer(a_imp_df$trttype),
  X = as.matrix(a_imp_df[, c_covs])
)

write_rds(c_trial_stan, "processed_data/trialsub_stan_data.rds")
