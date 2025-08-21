
source("scripts/00_config.R")
library(tidyverse)
library(multinma)
set.seed(123)

# Import pooled hazard ratios
a_imp_hr <- readRDS("outputs/conduct_model/cox_ma/ma_res.rds")

# Import synthetic IPD
a_imp_ipd <- readRDS("data/synth_ipd.rds")

# Format variables
a_ipd_fmt <- a_imp_ipd %>% 
  mutate(
    smoke = if_else(grepl("Former", smoke), 0L, 1L),
    race = as.factor(race),
    race = relevel(race, ref = "White"),
    across(c(comorbs, n_atc), as.numeric)
  )

# Get multivariate correlation
a_ipd_terms <- paste0("~ ", paste0(names(a_ipd_fmt[c(5:17)]), collapse = " + "))
a_ipd_mat <- model.matrix(as.formula(a_ipd_terms), data = a_ipd_fmt)[, -1]
a_ipd_cor <- cor(a_ipd_mat)

# Inverse CDF for categorical and ordinal variables in add_integration
qcat <- function(p, probs, values = seq_along(probs) - 1L) {
  probs <- probs / sum(probs)
  cdf <- cumsum(probs)
  p <- pmin(p, 1 - .Machine$double.eps)
  idx <- findInterval(p, c(0, cdf))
  values[idx]
}

# Dummy probabilities for ordinal and categorical variable sampling
comorb_prop <- c(
  "0" = 0.28, "1" = 0.41, "2" = 0.19, "3" = 0.07, "4" = 0.03, "5" = 0.01
)

conmed_prop <- c(
  "0" = 0.23, "1" = 0.25, "2" = 0.17, "3" = 0.07, "4" = 0.05, "5" = 0.05, 
  "6" = 0.04, "7" = 0.03, "8" = 0.03, "9" = 0.02, "10" = 0.02, "11" = 0.01, 
  "12" = 0.01, "13" = 0.01, "14" = 0.01
)

smoke_prop <- c("Former/Never" = 0.86, "Current" = 0.14)

race_prop <- c("Asian" = 0.3, "Black/AA" = 0.05, "Other" = 0.12, "White" = 0.54)

# Setup distribution parameters
b_params <- list(
  age10   = list(mu = 6, sd = 1),
  bmi     = list(mu = 5.4, sd = 1),
  hdl     = list(mu = 1.2, sd = 0.5),
  total   = list(mu = 3.8, sd = 1),
  durn    = list(shape = (1.4/1)^2, scale = (1^2)/1.4),
  egfr    = list(mu = 20, sd = 4),
  hba1c   = list(mu = 8, sd = 1),
  sbp     = list(mu = 40, sd = 3),
  sex     = list(prob = 0.54),
  smoke   = list(prob = 0.14),
  race    = list(probs = race_prop, values = seq_along(race_prop)),
  comorbs = list(probs = comorb_prop, values = as.integer(names(comorb_prop))),
  n_atc   = list(probs = conmed_prop, values = as.integer(names(conmed_prop)))
)

# Get sobol sequences for synthetic participants
b_sobol <- add_integration(
  data.frame(dummy = 1),
  age10 = distr(qnorm, b_params$age10$mu, b_params$age10$sd),
  bmi = distr(qnorm, b_params$bmi$mu, b_params$bmi$sd),
  hdl = distr(qnorm, b_params$hdl$mu, b_params$hdl$sd),
  total = distr(qnorm, b_params$total$mu, b_params$total$sd),
  comorbs = distr(qcat, probs = unname(b_params$comorbs$probs), values = b_params$comorbs$values),
  durn    = distr(qgamma, b_params$durn$shape, b_params$durn$scale),
  egfr    = distr(qnorm, b_params$egfr$mu, b_params$egfr$sd),
  hba1c   = distr(qnorm, b_params$hba1c$mu, b_params$hba1c$sd),
  n_atc   = distr(qcat, probs = unname(b_params$n_atc$probs), values = b_params$n_atc$values),
  raceAsian = distr(qbern, race_prop[[1]]),
  raceBlackAA = distr(qbern, race_prop[[2]]),
  raceOther = distr(qbern, race_prop[[3]]),
  sbp     = distr(qnorm, b_params$sbp$mu, b_params$sbp$sd),
  sex     = distr(qbern, b_params$sex$prob),
  smoke   = distr(qbern, b_params$smoke$prob),
  cor = a_ipd_cor,
  n_int = 64
)

b_sobol <- add_integration(
  b_df,
  age10 = distr(qnorm, age10mu, age10d),
  bmi = distr(qnorm, bmimu, bmid),
  hdl = distr(qnorm, hdlmu, hdld),
  total = distr(qnorm, totalmu, totald),
  comorbs = distr(qcat, probs = as.numeric(comorb_prop), values = as.integer(names(comorb_prop))),
  durn = distr(qgamma, durna, durnb),
  egfr = distr(qnorm, egfrmu, egfrd),
  hba1c = distr(qnorm, hba1cmu, hba1cd),
  n_atc = distr(qcat, probs = as.numeric(conmed_prop), values = as.integer(names(conmed_prop))),
  raceAsian = distr(qbern, asianp),
  raceBlackAA = distr(qbern, blackp),
  raceOther = distr(qbern, otherp),
  sbp = distr(qnorm, sbpmu, sbpd),
  sex = distr(qbern, sexp),
  smoke = distr(qbern, smokep),
  cor = a_ipd_cor,
  n_int = 64
)

# Extract sobol sequences for each covariate
b_sobol_ext <- b_sobol %>% 
  select(starts_with(".int")) %>% 
  rename_with(., ~ gsub(".int_", "", .)) %>% 
  unnest(cols = everything())

# Prepare pooled hazard ratios
b_hr <- a_imp_hr %>% 
  select(term, estimate) %>% 
  pivot_wider(names_from = "term", values_from = "estimate") %>% 
  rename_with(., ~ paste0(., "_hr"))

# Join hazard ratios with sobol sequences
b_join_hr <- b_sobol_ext %>% 
  cbind(b_hr) %>% 
  as_tibble() %>% 
  mutate(
    lp = age10*age10_hr + bmi*bmi_hr + hdl*chol_hdl_hr + total*chol_total_hr + 
      durn*durn_hr + egfr*egfr_hr + hba1c*hba1c_hr + sbp*sbp_hr + sex*sex_hr + 
      smoke*smoke_hr + raceAsian*raceasian_hr + raceBlackAA*raceblack_aa_hr + 
      raceOther * raceother_hr + comorbs*comorbs_hr + n_atc*n_atc_hr,
    attr = 1 / exp(lp)
  ) %>% 
  select(lp, attr, everything())

mean(b_join_hr$attr)

