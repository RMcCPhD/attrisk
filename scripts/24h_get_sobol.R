
# Obtain Sobol sequences
# Prepare trial-level matrices for CITL
source("scripts/00_config.R")
library(tidyverse)
library(multinma)
library(posterior)

# Import design features
a_imp_df <- readRDS("processed_data/latest_dataset.rds")

# Import marginals for baseline characteristics
c_imp_bl <- readRDS("processed_data/agg_bl_dists_ready.rds") %>% 
  select(trial_id, vars, values_scaled) %>% 
  pivot_wider(names_from = vars, values_from = values_scaled) %>% 
  mutate(male = 1 - male) %>% 
  rename(sex = male) %>% 
  group_by(trial_id) %>% 
  nest()

# Import correlation matrix
c_imp_cor <- readRDS("processed_data/citl_cor_mat.rds")

# Generate sobol sequences
c_sobol <- c_imp_bl %>% 
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
        black_aa    = distr(qbern, df$black),
        other       = distr(qbern, df$other),
        sbp         = distr(qnorm, df$sbp_mu, df$sbp_d),
        sex         = distr(qbern, df$sex),
        comorbs     = distr(qnbinom, size = df$comorbs_size, mu = df$comorbs_mu),
        n_atc       = distr(qnbinom, size = df$conmeds_size, mu = df$conmeds_mu),
        cor = c_imp_cor,
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
c_imp_hr <- readRDS("outputs/conduct_model/cox_ma/ma_res.rds") %>% 
  select(term, estimate) %>% 
  pivot_wider(names_from = "term", values_from = "estimate") %>% 
  rename_with(., ~ paste0(., "_loghr"))

# Rename betas
c_beta <- c(
  age10     = c_imp_hr$age10_loghr,
  bmi       = c_imp_hr$bmi_loghr,
  hdl       = c_imp_hr$chol_hdl_loghr,
  total     = c_imp_hr$chol_total_loghr,
  current   = c_imp_hr$smoke_loghr,
  t2dur     = c_imp_hr$durn_loghr,
  egfr      = c_imp_hr$egfr_loghr,
  hba1c     = c_imp_hr$hba1c_loghr,
  asian     = c_imp_hr$raceasian_loghr,
  black_aa  = c_imp_hr$raceblack_aa_loghr,
  other     = c_imp_hr$raceother_loghr,
  sbp       = c_imp_hr$sbp_loghr,
  sex       = c_imp_hr$sex_loghr,
  comorbs   = c_imp_hr$comorbs_loghr,
  n_atc     = c_imp_hr$n_atc_loghr
)

# Prepare CITL inputs ----------------------------------------------------------

# Parameters
I <- nrow(a_imp_df)  # N trials
J <- 128L            # N synthetic participants
Pnm <- names(c_beta) # Predictor names
pK <- length(Pnm)    # N predictors

# Construct large matrix of sobol sequences per trial (IxJ) x pK
Z_blocks <- vector("list", I)
for (i in seq_len(I)) {
  trial <- a_imp_df$trial_id[i]
  idx <- match(trial, c_sobol_tidy$trial_id)
  M <- c_sobol_tidy$sobol_ext[[idx]]
  Z_blocks[[i]] <- as.matrix(M)
}
Z <- do.call(rbind, Z_blocks)

# List of outputs
sobol_lst <- list(
  I = I,
  J = J,
  Pnm = Pnm,
  pK = pK,
  Z = Z,
  loghr = c_beta
)

str(sobol_lst)

# Save
saveRDS(sobol_lst, "processed_data/sobol_lst.rds")
