
rm(list = ls())
gc()

library(tidyverse)
library(pmsampsize)

### Aggregated models ----------------------------------------------------------

# Data
a_imp_stan <- readRDS("processed_data/design_stan_data.rds")
a_imp_agg <- read_csv("processed_data/data_both_models.csv")

# 34 covariates
glimpse(a_imp_agg)

# Sample size (number of arms across trials)
a_imp_agg %>% 
  group_by(trial_id) %>% 
  reframe(n = n()) %>% 
  count(n)

# Attrition prevalence
a_imp_agg %>%
  group_by(trial_id) %>% 
  reframe(
    n = sum(arm_n),
    attr = sum(arm_n_attr),
    attr_rate = attr/n
  ) %>% 
  reframe(mean_attr = mean(attr_rate))

# Test
pmsampsize(type = "b", csrsquared = 0.1, parameters = 14, prevalence = 0.152)

# Sample of 361 arms
# Define input grid
a_param_grid <- expand.grid(
  rsquared = seq(0.1, 0.2, by = 0.01),
  parameters = 1:34,
  prevalence = seq(0.1, 0.3, by = 0.05)
)

# Apply psampsize to each scenario and get outputs
a_sample <- a_param_grid %>% 
  mutate(
    sample = pmap(
      list(rsquared, parameters, prevalence),
      ~ pmsampsize(type = "b", csrsquared = ..1, parameters = ..2, prevalence = ..3)
    ),
    sample_size = map(sample, ~ .x$sample_size),
    events = map(sample, ~ .x$events),
    epp = map(sample, ~ .x$EPP),
    across(sample_size:epp, ~ unlist(.))
  )

# Find sample size closest to actual sample size
a_extract <- a_sample %>% 
  select(rsquared:prevalence, sample_size:epp) %>% 
  mutate(
    actual_sample = 498,
    actual_events = sum(a_imp_agg$arm_n_attr),
    sample_diff = abs(sample_size - actual_sample)
  ) %>% 
  arrange(sample_diff)

### IPD model ------------------------------------------------------------------

# Based on 41 T2DM trials in previous analysis
# 541 person-years of follow-up
# Mean follow-up of 1 year
# 14 predictors
# Event prevalence of 0.11

# Test
pmsampsize(
  type = "s", csrsquared = 0.1, parameters = 13, 
  timepoint = 1, rate = 0.11, meanfup = 1
)

# Sample of 361 arms
# Define input grid
b_param_grid <- expand.grid(
  rsquared = seq(0.1, 0.2, by = 0.01),
  parameters = 1:14,
  rate = 0.11,
  timepoint = 1,
  meanfup = 1
)

# Apply psampsize to each scenario and get outputs
b_sample <- b_param_grid %>% 
  mutate(
    sample = pmap(
      list(rsquared, parameters, rate, timepoint, meanfup),
      ~ pmsampsize(
        type = "s", csrsquared = ..1, parameters = ..2, rate = ..3,
        timepoint = ..4, meanfup = ..5
      )
    ),
    sample_size = map(sample, ~ .x$sample_size),
    events = map(sample, ~ .x$events),
    epp = map(sample, ~ .x$EPP),
    across(sample_size:epp, ~ unlist(.))
  )

# Find sample size closest to actual sample size
b_extract <- b_sample %>% 
  select(rsquared:rate, sample_size:epp) %>% 
  mutate(
    actual_sample = 541,
    actual_rate = 0.11,
    sample_diff = abs(sample_size - actual_sample)
  ) %>% 
  arrange(sample_diff)
