
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

# Drop a dummy variable in each group acting as a reference
# Dropping most common
a_imp_df %>% 
  select(n_250:reg_em) %>% 
  map_df(~ table(factor(.x, levels = c(0, 1))), .id = "variable") %>% 
  print(n = 34)

b_drop_ref <- a_imp_df %>% 
  select(-c(n_250, arm_100, dur_6m, arm_2, p3, bdbl, f_100, pom_1, reg_am))

# Design model dataset
names(b_drop_ref)
c_covs <- names(b_drop_ref[c(15:39)])

c_design_stan <- list(
  I = nrow(b_drop_ref),
  K = length(b_drop_ref),
  r = b_drop_ref$arm_n_attr,
  N = b_drop_ref$arm_n,
  dt = b_drop_ref$duration,
  arm = as.integer(b_drop_ref$trttype),
  X = as.matrix(b_drop_ref[, c_covs])
)

glimpse(c_design_stan)

# Save
write_rds(c_design_stan, "processed_data/design_stan_data.rds")

# Trial sub-model dataset
names(b_drop_ref)
d_covs <- names(b_drop_ref[c(10:39)])

d_trial_stan <- list(
  I = nrow(b_drop_ref),
  K = length(c_covs),
  r = b_drop_ref$arm_n_attr,
  N = b_drop_ref$arm_n,
  dt = b_drop_ref$duration,
  arm = as.integer(b_drop_ref$trttype),
  X = as.matrix(b_drop_ref[, d_covs])
)

write_rds(d_trial_stan, "processed_data/trialsub_stan_data.rds")

write_csv(
  a_imp_df %>% 
    select(arm_n_attr, arm_n, duration, trttype, d_covs) %>% 
    mutate(
      I = 361,
      K = 39
    ), 
  "processed_data/trialsub_stan_data.csv"
)
