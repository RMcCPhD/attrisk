
library(tidyverse)
library(pmsampsize)

## Aggregate-level models ######################################################

# Data
a_imp_stan <- readRDS("~/Documents/PhD/GitHub/attrisk/processed_data/stan_dat.rds")
a_imp_agg <- readRDS("~/Documents/PhD/GitHub/attrisk/processed_data/df_fu_n348.rds")

# 34 covariates (38 when including aggregated participant characteristics)
glimpse(a_imp_agg)

# Sample size (number of arms across trials)
a_imp_agg %>% 
  group_by(trial_id) %>% 
  reframe(n = n()) %>% 
  count(n)

# Sample of 369 arms
329 + (17*2) + (2*3)

# Define input grid
a_param_grid <- expand.grid(
  rsquared = seq(0.05, 0.2, by = 0.01),
  parameters = 1:38,
  prevalence = seq(0.1, 0.3, by = 0.05)
)

# Apply psampsize to each scenario
a_sample <- a_param_grid %>% 
  mutate(
    sample = pmap(
      list(rsquared, parameters, prevalence),
      ~ pmsampsize(type = "b", csrsquared = ..1, parameters = ..2, prevalence = ..3)
    )
  )

# Apply psampsize in each scenario
a_sample_estimates <- a_param_grid %>% 
  rowwise() %>% 
  mutate(
    sample_size = list(
      tryCatch({
        pmsampsize(
          type = "b",
          csrsquared = rsquared,
          parameters = parameters,
          prevalence = prevalence
        )
      },
      error = function(e) NA)
    )
  ) %>% 
  unnest(cols = c(sample_size))

