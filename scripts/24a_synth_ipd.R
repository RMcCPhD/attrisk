
source("scripts/00_config.R")
library(tidyverse)
library(multinma)
set.seed(123)

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

# Create synthetic IPD attrition data
a_ipd <- tibble(
  event = sample(rbinom(28, 1, prob = 0.2)),
  age10 = rnorm(28, mean = 6, sd = 1),
  bmi = rnorm(28, mean = 5.4, sd = 1),
  chol_hdl = rnorm(28, mean = 1.2, sd = 0.5),
  chol_total = rnorm(28, mean = 3.8, sd = 1),
  comorbs = sample(rep(names(comorb_prop), rmultinom(1, 28, prob = comorb_prop))),
  durn = rgamma(28, shape = (1.4/1)^2, scale = (1^2)/1.4),
  egfr = rnorm(28, mean = 20, sd = 4),
  hba1c = rnorm(28, mean = 8, sd = 1),
  n_atc = sample(rep(names(conmed_prop), rmultinom(1, 28, prob = conmed_prop))),
  race = sample(rep(names(race_prop), rmultinom(1, 28, prob = race_prop))),
  sbp = rnorm(28, mean = 40, sd = 3),
  sex = rbinom(28, 1, prob = 0.54),
  smoke = sample(rep(names(smoke_prop), rmultinom(1, 28, prob = smoke_prop)))
) %>% 
  mutate(
    id = row_number(),
    trial_id = ceiling(row_number() / 7),
    trttype = rbinom(28, 1, prob = 0.48)
  ) %>% 
  select(trial_id, id, trttype, event, everything())

# Save
saveRDS(a_ipd, "data/synth_ipd.rds")
