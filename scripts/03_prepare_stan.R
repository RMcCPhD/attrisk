
# Generate STAN-ready data for fitting design-stage model
source("scripts/00_config.R")
source("scripts/00_functions_packages_etc.R")

# Import tidied AgD data
imp_agd <- readRDS("processed_data/agd_steps_out.rds")$format_covars[[1]]

# Import exploratoryGompertz outputs
imp_expl <- readRDS("processed_data/shape_rate_design.rds")

# Import model matrix
imp_mat <- readRDS("processed_data/agd_covar_mat.rds")

# Create STAN data list
# I = Number of arms
# K = Number of covariates
# r = Arm-level attrition count
# N = Number of randomised arm participants
# dt = Follow-up duration in years
# lp = Shape linear predictor (shape ~ rate + sigma)
# beta_rate = Log-rate beta
# resid_sd = Residual SD
# X = Covariate matrix

stan_lst <- list(
  I = nrow(imp_agd),
  K = ncol(imp_mat),
  r = imp_agd$arm_attr,
  N = imp_agd$arm_n,
  dt = imp_agd$duration,
  beta_lograte = imp_expl$lograte_beta,
  alpha_shape = imp_expl$shape_alpha,
  resid_sd = imp_expl$uncertainty,
  X = imp_mat
)

str(stan_lst)
saveRDS(stan_lst, "processed_data/stan_data.rds")
