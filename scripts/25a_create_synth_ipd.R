
# Generate Sobol sequences
# 64 synthetic participants per aggregate trial
# Correlation matrix taken from nma_more_covariates

source("scripts/00_config.R")
source("scripts/00_paths.R")
library(tidyverse)
library(multinma)

# Import agg bl
a_imp_bl <- readRDS("processed_data/agg_bl_dists_ready.rds") %>% 
  select(trial_id, vars, values_scaled) %>% 
  pivot_wider(names_from = vars, values_from = values_scaled) %>% 
  mutate(male = 1 - male) %>% 
  rename(sex = male) %>% 
  group_by(trial_id) %>% 
  nest()

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
c_sobol <- a_imp_bl %>% 
  mutate(
    sobol = map(data, function(df) {
      add_integration(
        df,
        age10       = distr(qnorm, df$age_mu, df$age_d),
        bmi         = distr(qnorm, df$bmi_mu, df$bmi_d),
        hdl         = distr(qnorm, df$hdl_mu, df$hdl_d),
        total       = distr(qnorm, df$total_chol_mu, df$total_chol_d),
        current     = distr(qbern, df$current),
        t2dur       = distr(qgamma, df$t2dur_mu, df$t2dur_d),
        egfr        = distr(qnorm, df$egfr_mu, df$egfr_d),
        hba1c       = distr(qnorm, df$hba1c_mu, df$hba1c_d),
        asian       = distr(qbern, df$asian),
        black       = distr(qbern, df$black),
        other       = distr(qbern, df$other),
        sbp         = distr(qnorm, df$sbp_mu, df$sbp_d),
        sex         = distr(qbern, df$sex),
        # comorbs   = distr(qcat, probs = unname(atc$n_atc), values = atc$prop),
        # n_atc     = distr(qcat, probs = unname(atc$n_atc), values = atc$prop),
        cor = as.matrix(b_cor_rm),
        n_int = 64
      )
    })
  )

# Tidy sequences
c_sobol_tidy <- c_sobol %>% 
  mutate(
    sobol_ext = map(sobol, function(df) {
      df %>% 
        select(starts_with(".int")) %>% 
        rename_with(., ~ gsub(".int_", "", .)) %>% 
        unnest(cols = everything())
    })
  )

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
      lp <- as.vector(X %*% d_beta)
      tibble(lp = lp, lp_hr = exp(lp), risk = 1 - exp(-0.153 * lp_hr))
    }),
    avg_risk = map_dbl(lp_tbl, ~ sum(.x$risk) / 64)
  ) %>% 
  select(trial_id, lp_tbl, avg_risk)
