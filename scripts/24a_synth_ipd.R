
source("scripts/00_config.R")
library(tidyverse)
library(multinma)
set.seed(123)

# List of distribution inputs, i.e. dummy probabilities and parameters
# Scaled as done for the Cox models
# Durn: Prepared for gamma distribution (i.e. shape, scale; right-skewed)
# Comorbs: 0 to 5
# ATC: 0 to 10
# Smoke: Former/Never, Current
# Race: Asian, Black/AA, Other, White
a_inputs <- list(
  age10   = c(6, 1),
  sex     = c(0.54),
  bmi     = c(5.4, 1),
  hdl     = c(1.2, 0.5),
  total   = c(3.8, 1),
  comorbs = c(0.28, 0.41, 0.19, 0.07, 0.03, 0.01),
  durn    = c((1.4/1)^2, (1^2)/1.4),
  egfr    = c(20, 4),
  hba1c   = c(8, 1),
  sbp     = c(40, 3),
  n_atc   = c(0.23, 0.25, 0.17, 0.07, 0.07, 0.05, 0.05, 0.04, 0.03, 0.03, 0.01),
  smoke   = c(0.86, 0.14),
  race    = c(0.3, 0.05, 0.12, 0.54)
)

vec_race <- c("asian", "black_aa", "other", "white")

# Create synthetic IPD groups - 64 patients, 4 trials
a_ipd_grps <- tibble(idx = 1:256) %>% 
  mutate(trial_id = ceiling(row_number() / 64)) %>% 
  group_by(trial_id) %>% 
  mutate(id = row_number()) %>% 
  ungroup() %>% 
  select(-idx) %>% 
  group_by(trial_id) %>% 
  nest() %>% 
  ungroup()

# Simulate variables for each patient per trial
b_synth <- a_ipd_grps %>% 
  mutate(
    data = map(data, function(df) {
      df %>% 
        mutate(
          age10   = rnorm(64, a_inputs$age10[[1]], a_inputs$age10[[2]]),
          sex     = rbinom(64, 1, a_inputs$sex),
          bmi     = rnorm(64, a_inputs$bmi[[1]], a_inputs$bmi[[2]]),
          hdl     = rnorm(64, a_inputs$hdl[[1]], a_inputs$hdl[[2]]),
          total   = rnorm(64, a_inputs$total[[1]], a_inputs$total[[2]]),
          comorbs = sample(rep(0:5, rmultinom(1, 64, a_inputs$comorbs))),
          durn    = rgamma(64, shape = a_inputs$durn[[1]], scale = a_inputs$durn[[2]]),
          egfr    = rnorm(64, a_inputs$egfr[[1]], a_inputs$egfr[[2]]),
          hba1c   = rnorm(64, a_inputs$hba1c[[1]], a_inputs$hba1c[[2]]),
          sbp     = rnorm(64, a_inputs$sbp[[1]], a_inputs$sbp[[2]]),
          n_atc   = sample(rep(0:10, rmultinom(1, 64, a_inputs$n_atc))),
          smoke   = sample(rep(0:1, rmultinom(1, 64, a_inputs$smoke))),
          race    = sample(rep(vec_race, rmultinom(1, 64, a_inputs$race)))
        )
    })
  ) %>% 
  unnest(data)

# Summarise distributions per trial - numeric
b_synth_sum <- b_synth %>% 
  group_by(trial_id) %>% 
  reframe(
    across(
      c(age10:total, durn:sbp), 
      ~ unlist(list(mean = mean(.), sd = sd(.)))
    )
  ) %>% 
  cbind(tibble(marg = c("mean", "sd"))) %>% 
  as_tibble()

# Summarise categorical variables
b_synth_cat <- list(
  comorbs = b_synth %>% group_by(trial_id) %>% count(comorbs) %>% mutate(prop = n/64),
  n_atc   = b_synth %>% group_by(trial_id) %>% count(n_atc) %>% mutate(prop = n/64),
  race    = b_synth %>% group_by(trial_id) %>% count(race) %>% mutate(prop = n/64),
  smoke   = b_synth %>% group_by(trial_id) %>% count(smoke) %>% mutate(prop = n/64),
  sex     = b_synth %>% group_by(trial_id) %>% count(sex) %>% mutate(prop = n/64)
)

# Save
saveRDS(b_synth, "data/synth_ipd.rds")
saveRDS(b_synth_sum, "data/synth_ipd_sum.rds")
saveRDS(b_synth_cat, "data/synth_ipd_sum_cat.rds")
