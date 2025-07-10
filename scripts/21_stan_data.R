
source("scripts/00_config.R")
library(tidyverse)

# Import data and model coefficients
a_imp_df <- readRDS("processed_data/latest_dataset.rds")
a_imp_mdl_design_nm <- read_csv("outputs/model_shape/design/lm_res.csv")
a_imp_mdl_design <- readRDS("outputs/model_shape/design/lm_res_shrinkage.rds")
a_imp_mdl_conduct_nm <- read_csv("outputs/model_shape/conduct/lm_res.csv")
a_imp_mdl_conduct <- readRDS("outputs/model_shape/conduct/lm_res_shrinkage.rds")

# Design data ------------------------------------------------------------------

# Remove baseline characteristics
b_design_df <- a_imp_df %>% 
  select(-c(contains("dbl"), bmi_lower:bmi_upper, male_prop:egfr_mean))

# Create model matrix for predictors
b_vars <- names(b_design_df[9:25])
b_fml <- as.formula(paste0("~ ", paste0(b_vars, collapse = " + ")))
b_design_mat <- model.matrix(b_fml, data = b_design_df)[, -1]

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
  I = nrow(b_design_df),
  K = ncol(b_design_mat),
  r = b_design_df$attr,
  N = b_design_df$arm_n,
  dt = b_design_df$duration,
  B = b_coefs_mat,
  X = b_design_mat
)

# Save
saveRDS(b_design_stan, "processed_data/latest_stan_design.rds")

# Conduct data -----------------------------------------------------------------

# Same dataset, bl characteristics not included as predictors here
c_conduct_df <- b_design_df

# Create model matrix for predictors
c_vars <- names(c_conduct_df[c(9:25)])
c_fml <- as.formula(paste0("~ ", paste0(c_vars, collapse = " + ")))
c_conduct_mat <- model.matrix(c_fml, data = c_conduct_df)[, -1]

# Prepare conduct coefficients
c_conduct_coefs <- a_imp_mdl_conduct_nm %>% 
  filter(!(term %in% c("(Intercept)", "rate"))) %>% 
  select(term, estimate) %>% 
  pivot_wider(names_from = "term", values_from = "estimate")

c_coefs_names <- names(c_conduct_coefs)
c_coefs_fml <- as.formula(paste0("~ ", paste0(c_coefs_names, collapse = " + ")))

c_coefs_beta <- tibble(term = c_coefs_names, estimate = a_imp_mdl_conduct) %>% 
  pivot_wider(names_from = "term", values_from = "estimate")

c_coefs_mat <- model.matrix(c_coefs_fml, data = c_coefs_beta)[-1]

# Prepare STAN data
c_conduct_stan <- list(
  I = nrow(c_conduct_df),
  K = ncol(c_conduct_mat),
  r = c_conduct_df$attr,
  N = c_conduct_df$arm_n,
  dt = c_conduct_df$duration,
  B = c_coefs_mat,
  X = c_conduct_mat
)

# Save
saveRDS(c_conduct_stan, "processed_data/latest_stan_conduct.rds")
