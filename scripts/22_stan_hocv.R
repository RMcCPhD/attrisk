
source("scripts/00_config.R")
library(tidyverse)

# Import data and model coefficients
a_imp_df <- readRDS("processed_data/latest_dataset.rds")
a_imp_mdl_design_nm <- read_csv("outputs/model_shape/design/lm_res.csv")
a_imp_mdl_design <- readRDS("outputs/model_shape/design/lm_res_shrinkage.rds")

# Design data ------------------------------------------------------------------

# Remove baseline characteristics
b_design_df <- a_imp_df %>% 
  select(-c(contains("dbl"), bmi_lower:bmi_upper, male_prop:egfr_mean))

# Trial IDs where there are more than one arm
b_mt1_arm <- b_design_df %>% 
  group_by(trial_id) %>% 
  filter(length(unique(trttype)) > 1) %>% 
  ungroup() %>% 
  distinct(trial_id)

mt1_trials <- b_mt1_arm$trial_id

# Remove 14 trials with >1 treatment arm
b_rm14 <- b_design_df %>% filter(!trial_id %in% mt1_trials)

# Create model matrix for predictors
b_vars <- names(b_rm14[9:25])
b_fml <- as.formula(paste0("~ ", paste0(b_vars, collapse = " + ")))
b_design_mat <- model.matrix(b_fml, data = b_rm14)[, -1]

# Prepare design coefficients
b_design_coefs <- a_imp_mdl_design_nm %>% 
  filter(!(term %in% c("(Intercept)", "rate"))) %>% 
  select(term, estimate) %>% 
  pivot_wider(names_from = "term", values_from = "estimate")

b_coefs_names <- names(b_design_coefs)
b_coefs_fml <- as.formula(paste0("~ ", paste0(b_coefs_names, collapse = " + ")))

b_coefs_beta <- tibble(term = b_coefs_names, estimate = a_imp_mdl_design) %>% 
  pivot_wider(names_from = "term", values_from = "estimate")

b_coefs_mat <- model.matrix(b_coefs_fml, data = b_coefs_beta)[-1]

# Prepare STAN data
b_design_stan <- list(
  I = nrow(b_rm14),
  K = ncol(b_design_mat),
  r = b_rm14$attr,
  N = b_rm14$arm_n,
  dt = b_rm14$duration,
  B = b_coefs_mat,
  X = b_design_mat
)

str(b_design_stan)

# Save
saveRDS(b_design_stan, "processed_data/design_hocv.rds")
