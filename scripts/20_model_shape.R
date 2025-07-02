
source("scripts/00_config.R")
library(tidyverse)

# Import dataset
a_imp_df <- readRDS("processed_data/latest_dataset.rds")

# Trials of interest
z_ipd_trials <- a_imp_df %>% filter(source == "ipd") %>% distinct(trial_id)
z_ipd_trials <- as.vector(z_ipd_trials$trial_id)

# Import vivli intercept-only model results
a_imp_mdl <- read_csv("data/attrition_intercept_models/tidy_res.csv")

# Import intercept-only models from attrition rates paper
a_imp_mdl_rates <- read_csv("data/attrition_rates_data/parameters.csv")

# Prepare newer model results
# Was modelled in years, no scaling needed
b_mdl_res <- a_imp_mdl %>% 
  rename(
    trial_id = nct_id,
    shape = shape_est,
    rate = rate_est
  ) %>% 
  select(-c(shape_se, rate_se))

# All present as expected
sum(b_mdl_res$trial_id %in% z_ipd_trials)

# Prepare rates model data
# Remove trials in newer set
# Scale to duration in years
c_tidy_rates <- a_imp_mdl_rates %>% 
  filter(dist == "gompertz" & parameter %in% c("rate", "shape")) %>% 
  select(-c(dist, se)) %>% 
  pivot_wider(names_from = "parameter", values_from = "est") %>% 
  rename(trial_id = ctgov) %>% 
  filter(
    !(trial_id %in% b_mdl_res$trial_id),
    trial_id %in% z_ipd_trials
  ) %>% 
  mutate(
    shape = shape * 365,
    rate = rate + log(365)
  )

# 4 additional trials present
sum(c_tidy_rates$trial_id %in% z_ipd_trials)

# Join model datasets, 76 trials
d_join_mdl <- b_mdl_res %>% full_join(c_tidy_rates)

# Join characteristics
e_join_chars <- d_join_mdl %>% left_join(a_imp_df)

# Design model -----------------------------------------------------------------

# Data
f_design_df <- e_join_chars %>% select(-contains("mean"))

# Prepare model formula
z_covars <- paste0(names(f_design_df[8:34]), collapse = " + ")
z_fml <- as.formula(paste0("shape ~ rate + ", z_covars))

# Correlation beteween rate and shape
# -0.453
cor(f_design_df$shape, f_design_df$rate, method = "pearson")

# Fit linear model
f_mdl <- lm(z_fml, data = f_design_df)
summary(f_mdl)

# Conduct model ----------------------------------------------------------------
