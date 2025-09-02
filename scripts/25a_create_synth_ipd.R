
# Generate Sobol sequences
# 64 synthetic participants per aggregate trial
# Correlation matrix taken from nma_more_covariates

source("scripts/00_config.R")
source("scripts/00_paths.R")
library(tidyverse)
library(multinma)

# Import agg bl
a_imp_bl <- readRDS("processed_data/agg_bl_dists_ready.rds") %>% 
  select(vars, values_scaled) %>% 
  group_by(vars) %>% 
  nest() %>% 
  mutate(data = map_dbl(data, ~ .x$values_scaled)) %>% 
  pivot_wider(names_from = vars, values_from = data) %>% 
  mutate(male = 1 - male) %>% 
  rename(sex = male)

# Prepare correlation matrix ---------------------------------------------------

# Import covariate correlation
b_imp_cor <- readRDS(z_imp_cor)

# Row and column names
b_cor_nms <- c(
  "age10", "rm", "bmi", "rm", "hdl", "rm", "total_chol", "current", "rm",
  "t2dur", "egfr", "hba1c", "rm", "rm", "rm", "rm", "rm", "asian", "black",
  "other", "rm", "sbp", "sex", "rm", "rm"
)

rownames(b_imp_cor) <- b_cor_nms
colnames(b_imp_cor) <- b_cor_nms

# Remove external variables
b_cor_rm <- b_imp_cor[!rownames(b_imp_cor) == "rm", !rownames(b_imp_cor) == "rm"]

# Generate Sobol sequences -----------------------------------------------------

# Helper: Inverse CDF for con meds and comorbs
qcat <- function(p, probs, values = seq_along(probs) - 1L) {
  probs <- probs / sum(probs)
  cdf <- cumsum(probs)
  p <- pmin(p, 1 - .Machine$double.eps)
  idx <- findInterval(p, c(0, cdf))
  values[idx]
}

# Generate sobol sequences
c_sobol <- add_integration(
  a_imp_bl,
  age10       = distr(qnorm, a_imp_bl$age_mu, a_imp_bl$age_d),
  bmi         = distr(qnorm, a_imp_bl$bmi_mu, a_imp_bl$bmi_d),
  hdl         = distr(qnorm, a_imp_bl$hdl_mu, a_imp_bl$hdl_d),
  total       = distr(qnorm, a_imp_bl$total_chol_mu, a_imp_bl$total_chol_d),
  current     = distr(qbern, a_imp_bl$current),
  t2dur       = distr(qgamma, a_imp_bl$t2dur_mu, a_imp_bl$t2dur_d),
  egfr        = distr(qnorm, a_imp_bl$egfr_mu, a_imp_bl$egfr_d),
  hba1c       = distr(qnorm, a_imp_bl$hba1c_mu, a_imp_bl$hba1c_d),
  asian       = distr(qbern, a_imp_bl$asian),
  black       = distr(qbern, a_imp_bl$black),
  other       = distr(qbern, a_imp_bl$other),
  sbp         = distr(qnorm, a_imp_bl$sbp_mu, a_imp_bl$sbp_d),
  sex         = distr(qbern, a_imp_bl$sex),
  # comorbs   = distr(qcat, probs = unname(atc$n_atc), values = atc$prop),
  # n_atc       = distr(qcat, probs = unname(atc$n_atc), values = atc$prop),
  cor = b_cor_rm,
  n_int = 64
)

# Tidy sequences
c_sobol_tidy <- c_sobol %>% 
  select(starts_with(".int")) %>% 
  rename_with(., ~ gsub(".int_", "", .)) %>% 
  unnest(cols = everything())

# Import log-hazard estimates
d_imp_hr <- readRDS("outputs/conduct_model/cox_ma/ma_res.rds") %>% 
  select(term, estimate) %>% 
  pivot_wider(names_from = "term", values_from = "estimate") %>% 
  rename_with(., ~ paste0(., "_loghr"))

# Rename betas
d_beta <- c(
  age10     = d_imp_hr$age10_loghr,
  bmi       = d_imp_hr$bmi_loghr,
  hdl       = d_imp_hr$chol_hdl_loghr,
  total     = d_imp_hr$chol_total_loghr,
  current   = d_imp_hr$smoke_loghr,
  t2dur     = d_imp_hr$durn_loghr,
  egfr      = d_imp_hr$egfr_loghr,
  hba1c     = d_imp_hr$hba1c_loghr,
  asian     = d_imp_hr$raceasian_loghr,
  black_aa  = d_imp_hr$raceblack_aa_loghr,
  other     = d_imp_hr$raceother_loghr,
  sbp       = d_imp_hr$sbp_loghr,
  sex       = d_imp_hr$sex_loghr
  # comorbs   = d_imp_hr$comorbs_loghr,
  # n_atc     = d_imp_hr$n_atc_loghr,
)

# Compute linear predictor per trial
# Get average predicted attrition risk per trial
d_lp <- c_sobol_tidy %>% 
  mutate(
    lp_tbl = map(sobol_ext, function(df) {
      X <- as.matrix(df)
      lp <- as.vector(X %*% b_beta)
      tibble(lp = lp, lp_hr = exp(lp), risk = 1 - exp(-0.153 * lp_hr))
    }),
    avg_risk = map_dbl(lp_tbl, ~ sum(.x$risk) / 64)
  ) %>% 
  select(trial_id, lp_tbl, avg_risk)