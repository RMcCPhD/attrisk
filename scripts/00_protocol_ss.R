
# Approximate sample size adequacy for AgD and IPD models

rm(list = ls())
gc()

library(tidyverse)
library(pmsampsize)

# Helper functions for running psampsize (AgD: Binary, IPD: Survival)
pms_binary <- function(csrsquared, parameters, prevalence) {
  
  out <- pmsampsize(
    type = "b",
    csrsquared = csrsquared,
    parameters = parameters,
    prevalence = prevalence
  )
  
  tibble(
    csrsquared = csrsquared,
    parameters = parameters,
    prevalence = prevalence,
    sample_size = out$sample_size,
    events = out$events,
    epp = out$EPP
  )
  
}

pms_surv <- function(csrsquared, parameters, rate, timepoint, meanfup) {
  
  out <- pmsampsize(
    type = "s",
    csrsquared = csrsquared,
    parameters = parameters,
    rate = rate,
    timepoint = timepoint,
    meanfup = meanfup
  )
  
  tibble(
    csrsquared = csrsquared,
    parameters = parameters,
    rate = rate,
    timepoint = timepoint,
    meanfup = meanfup,
    sample_size = out$sample_size,
    events = out$events,
    epp = out$EPP
  )
  
}

# AgD model --------------------------------------------------------------------

# Data
a_imp_stan <- readRDS("processed_data/design_stan_data.rds")
a_imp_agd <- readRDS("processed_data/latest_dataset.rds")

# Define observed sample (344)
a_agd_n <- nrow(a_imp_agd)

# Attrition prevalence (median: 0.127, overall 0.169)
a_agd_attr_fu <- a_imp_agd %>% 
  select(trial_id, arm_n, attr, duration) %>% 
  mutate(rate = attr / arm_n) %>% 
  reframe(
    med_rate = median(rate),
    med_fu = mean(duration)
  )

a_agd_attr_overall <- sum(a_imp_agd$attr) / sum(a_imp_agd$arm_n)

# 25 covariates (max number of covariates including dummy coded factors)
glimpse(a_imp_stan)

# Create expansion grid to explore scenarios
a_agd_grid <- expand_grid(
  csrsquared = seq(0.05, 0.20, by = 0.01),
  parameters = 1:25,
  prevalence = c(
    round(a_agd_attr_overall, 3),
    0.10, 0.15, 0.20, 0.25, 0.30
  ) %>% unique() %>% sort()
)

# Run pmsampsize and extract results
a_agd_res <- pmap_dfr(
  a_agd_grid,
  \(csrsquared, parameters, prevalence) {
    pms_binary(
      csrsquared = csrsquared,
      parameters = parameters,
      prevalence = prevalence
    )
  }
)

a_agd_out <- a_agd_res %>%
  mutate(
    available_sample = a_agd_n,
    feasible = sample_size <= available_sample,
    sample_diff = available_sample - sample_size
  )

# Get largest feasible model for each scenario
a_agd_sum <- a_agd_out %>%
  filter(feasible) %>%
  group_by(csrsquared, prevalence) %>%
  slice_max(order_by = parameters, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(prevalence, csrsquared)

# Scenario closest to observed data
a_agd_closest <- a_agd_sum %>% 
  filter(
    abs(prevalence - a_agd_attr_overall) == min(abs(prevalence - a_agd_attr_overall))
  ) %>%
  arrange(csrsquared)

# Conservative scenarios
agd_consv <- a_agd_sum %>%
  summarise(
    min_feasible_parameters = min(parameters),
    max_feasible_parameters = max(parameters),
    min_epp = min(epp),
    max_epp = max(epp)
  )

# IPD model --------------------------------------------------------------------

# Data from previous attrition rates analysis
b_imp_ipd <- read_csv("../attrition_rates_modelling/usable_data/total_attrition/study_chars.csv") %>% 
  filter(condition == "Type 2 diabetes")

# Sample size (28335)
b_ipd_n <- sum(b_imp_ipd$n)

# Median event rate and follow-up (0.107, 1.10)
b_ipd_rate_fu <- b_imp_ipd %>% 
  mutate(
    rate = total_events / n,
    fu_yr = fu_mth / 12
  ) %>% 
  reframe(
    med_rate = median(rate),
    med_fu = mean(fu_yr)
  )

# Rounding to 1 median year of FU for simplicity
b_ipd_rate <- 0.107
b_ipd_fu <- 1
b_ipd_tmpt <- 1
b_ipd_covars <- 13

# Expansion grid to explore scenarios
b_ipd_grid <- expand_grid(
  csrsquared = seq(0.05, 0.20, by = 0.01),
  parameters = 1:13,
  rate = c(b_ipd_rate, 0.10, 0.15, 0.20) %>% unique() %>% sort(),
  timepoint = b_ipd_tmpt,
  meanfup = b_ipd_fu
)

# Run pmsampsize and extract results
b_ipd_res <- pmap_dfr(
  b_ipd_grid,
  \(csrsquared, parameters, rate, timepoint, meanfup) {
    pms_surv(
      csrsquared = csrsquared,
      parameters = parameters,
      rate = rate,
      timepoint = timepoint,
      meanfup = meanfup
    )
  }
)

b_ipd_out <- b_ipd_res %>% 
  mutate(
    available_sample = b_ipd_n,
    feasible = sample_size <= available_sample,
    sample_diff = available_sample - sample_size
  )

# Get largest feasible model
b_ipd_sum <- b_ipd_out %>% 
  filter(feasible) %>%
  group_by(csrsquared, rate, timepoint, meanfup) %>%
  slice_max(order_by = parameters, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(rate, csrsquared)

# Get scenario matching observed data
b_ipd_closest <- b_ipd_out %>% 
  filter(rate == b_ipd_rate) %>%
  arrange(csrsquared)

# Convervative scenarios
b_ipd_consv <- b_ipd_out %>% 
  summarise(
    min_feasible_parameters = min(parameters),
    max_feasible_parameters = max(parameters),
    min_epp = min(epp),
    max_epp = max(epp)
  )

# Get reportable scenario ------------------------------------------------------

c_target_csrq <- 0.10 # Conservative

c_agd <- a_agd_out %>% 
  filter(
    feasible,
    abs(prevalence - a_agd_attr_overall) == min(abs(prevalence - a_agd_attr_overall)),
    csrsquared == c_target_csrq
  ) %>%
  slice_max(order_by = parameters, n = 1, with_ties = FALSE)

c_ipd <- b_ipd_out %>%
  filter(
    feasible,
    rate == b_ipd_rate,
    csrsquared == c_target_csrq
  ) %>%
  slice_max(order_by = parameters, n = 1, with_ties = FALSE)

write_csv(b_ipd_res, "../../Attrisk/GitHub/pmsampsize_ipd_all.csv")
write_csv(c_ipd, "../../Attrisk/GitHub/pmsampsize_ipd_best.csv")
