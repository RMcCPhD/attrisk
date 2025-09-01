
source("scripts/00_config.R")
library(tidyverse)
library(multinma)
set.seed(123)

# Import synthetic IPD and summaries
a_imp_ipd <- readRDS("data/synth_ipd.rds")
a_imp_sum <- readRDS("data/synth_ipd_sum.rds")
a_imp_cat <- readRDS("data/synth_ipd_sum_cat.rds")

# Format race, nest data
a_ipd_fmt <- a_imp_ipd %>% 
  mutate(
    race = as.factor(race),
    race = relevel(race, ref = "white"),
  ) %>% 
  group_by(trial_id) %>% 
  nest() %>% 
  ungroup()

# Get correlation matrix for each trial
b_ipd_cor <- a_ipd_fmt %>% 
  mutate(
    cor = map(data, function(df) {
      terms <- paste0("~ ", paste0(names(df[c(2:14)]), collapse = " + "))
      fml <- model.matrix(as.formula(terms), data = df)[, -1]
      cor(fml)
    })
  )

# Add summaries to nested data
b_add_sum <- b_ipd_cor %>% 
  left_join(
    a_imp_sum %>% 
      pivot_wider(names_from = marg, values_from = age10:sbp) %>% 
      group_by(trial_id) %>%
      nest(.key = "numeric") %>% 
      ungroup()
  ) %>% 
  left_join(
    a_imp_cat %>% 
      bind_rows(.id = "var") %>% 
      group_by(trial_id, var) %>% 
      nest() %>% 
      ungroup() %>% 
      pivot_wider(names_from = var, values_from = data) %>% 
      rename(conmeds = n_atc)
  )

# Inverse CDF for categorical and ordinal variables in add_integration
qcat <- function(p, probs, values = seq_along(probs) - 1L) {
  probs <- probs / sum(probs)
  cdf <- cumsum(probs)
  p <- pmin(p, 1 - .Machine$double.eps)
  idx <- findInterval(p, c(0, cdf))
  values[idx]
}

# Generate sobol sequences
b_sobol <- b_add_sum %>% 
  mutate(
    sobol = pmap(
      list(cor, numeric, comorbs, conmeds, race, smoke, sex),
      function(r, df, morb, atc, rce, smk, male) {
        add_integration(
          df,
          age10       = distr(qnorm, df$age10_mean, df$age10_sd),
          sex         = distr(qbern, male$prop[[2]]),
          bmi         = distr(qnorm, df$bmi_mean, df$bmi_sd),
          hdl         = distr(qnorm, df$hdl_mean, df$hdl_sd),
          total       = distr(qnorm, df$total_mean, df$total_sd),
          comorbs     = distr(qcat, probs = unname(morb$comorbs), values = morb$prop),
          durn        = distr(qgamma, df$durn_mean, df$durn_sd),
          egfr        = distr(qnorm, df$egfr_mean, df$egfr_sd),
          hba1c       = distr(qnorm, df$hba1c_mean, df$hba1c_sd),
          sbp         = distr(qnorm, df$sbp_mean, df$sbp_sd),
          n_atc       = distr(qcat, probs = unname(atc$n_atc), values = atc$prop),
          smoke       = distr(qbern, smk$prop[[2]]),
          asian       = distr(qbern, rce$prop[[1]]),
          black_aa    = distr(qbern, rce$prop[[2]]),
          other       = distr(qbern, rce$prop[[3]]),
          cor = r,
          n_int = 64
        )
      }
    ),
    sobol_ext = map(sobol, function(df) {
      df %>% 
        select(starts_with(".int")) %>% 
        rename_with(., ~ gsub(".int_", "", .)) %>% 
        unnest(cols = everything())
    })
  )

# Import log-hazard estimates
b_imp_hr <- readRDS("outputs/conduct_model/cox_ma/ma_res.rds") %>% 
  select(term, estimate) %>% 
  pivot_wider(names_from = "term", values_from = "estimate") %>% 
  rename_with(., ~ paste0(., "_loghr"))

# Rename betas
b_beta <- c(
  age10     = b_imp_hr$age10_loghr,
  sex       = b_imp_hr$sex_loghr,
  bmi       = b_imp_hr$bmi_loghr,
  hdl       = b_imp_hr$chol_hdl_loghr,
  total     = b_imp_hr$chol_total_loghr,
  comorbs   = b_imp_hr$comorbs_loghr,
  durn      = b_imp_hr$durn_loghr,
  egfr      = b_imp_hr$egfr_loghr,
  hba1c     = b_imp_hr$hba1c_loghr,
  sbp       = b_imp_hr$sbp_loghr,
  n_atc     = b_imp_hr$n_atc_loghr,
  smoke     = b_imp_hr$smoke_loghr,
  asian     = b_imp_hr$raceasian_loghr,
  black_aa  = b_imp_hr$raceblack_aa_loghr,
  other     = b_imp_hr$raceother_loghr
)

# Compute linear predictor per trial
# Get average predicted attrition risk per trial
c_lp <- b_sobol %>% 
  mutate(
    lp_tbl = map(sobol_ext, function(df) {
      X <- as.matrix(df)
      lp <- as.vector(X %*% b_beta)
      tibble(lp = lp, lp_hr = exp(lp), risk = 1 - exp(-0.153 * lp_hr))
    }),
    avg_risk = map_dbl(lp_tbl, ~ sum(.x$risk) / 64)
  ) %>% 
  select(trial_id, lp_tbl, avg_risk)
